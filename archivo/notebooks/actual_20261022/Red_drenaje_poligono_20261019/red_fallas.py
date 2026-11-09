"""Red D8 completa, ksn móvil y perfiles que cruzan trazas cartográficas."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
from collections import deque
import json,hashlib
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.windows import Window
from scipy.spatial import cKDTree
from shapely.geometry import LineString,Point
import shapely

CRS='EPSG:32718'
OFF={1:(0,1),2:(1,1),4:(1,0),8:(1,-1),16:(0,-1),32:(-1,-1),64:(-1,0),128:(-1,1)}

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def leer_poligono(path):
    g=gpd.read_file(path)
    if g.crs is None:raise ValueError('Polígono sin CRS.')
    g=g.to_crs(CRS)
    if not g.is_valid.all():raise ValueError('Polígono inválido; revisar la geometría original.')
    if not g.geom_type.isin(['Polygon','MultiPolygon']).all():raise ValueError('La capa debe contener polígonos.')
    return g,g.geometry.union_all()

def extraer_grafo(dem,area,flow,amin,cache,theta=.4,A0=1e6,min_elevation=1e-6):
    """No se recalcula flujo en un DEM recortado: se conserva el aporte exterior al AOI."""
    cache=Path(cache);cache.mkdir(parents=True,exist_ok=True)
    signature={'algoritmo':'D8_sparse_v2_tierra','area_minima_m2':float(amin),'theta':theta,'A0':A0,
               'elevacion_minima_dem_m':min_elevation,
               'fuentes':{str(p):sha(p) for p in [dem,area,flow]}}
    key=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:18]
    target=cache/f'grafo_{key}.npz';meta=cache/f'grafo_{key}.json'
    if target.exists() and meta.exists():
        print('Recuperando grafo completo de caché',flush=True)
        with np.load(target,allow_pickle=False) as z:g={k:z[k] for k in z.files}
        info=json.loads(meta.read_text(encoding='utf8'));info['cache']=str(target);return g,info
    ids=[];areas=[];zz=[];fds=[]
    with rasterio.open(area) as a,rasterio.open(flow) as f,rasterio.open(dem) as z:
        if any(s.crs!=rasterio.crs.CRS.from_epsg(32718) for s in [a,f,z]):raise ValueError('CRS incompatible.')
        if any((s.transform!=a.transform or s.shape!=a.shape) for s in [f,z]):raise ValueError('DEM, D8 y área no alineados.')
        H,W=a.shape;tr=a.transform
        if tr.b or tr.d:raise ValueError('Raster rotado no soportado.')
        for top in range(0,H,256):
            win=Window(0,top,W,min(256,H-top))
            av=a.read(1,window=win);fv=f.read(1,window=win);zv=z.read(1,window=win,masked=True).astype(float).filled(np.nan)
            mask=np.isfinite(av)&(av>=amin)&np.isfinite(zv)
            if min_elevation is not None:mask &= zv>min_elevation
            rr,cc=np.nonzero(mask)
            ids.append(((rr+top)*W+cc).astype(np.int64));areas.append(av[mask]);zz.append(zv[mask]);fds.append(fv[mask])
            if top%2048==0:print(f'Lectura global D8/área: fila {top}/{H}',flush=True)
    gid=np.concatenate(ids);A=np.concatenate(areas);Z=np.concatenate(zz);fd=np.concatenate(fds);N=len(gid)
    rr=gid//W;cc=gid%W
    dr=np.zeros(256,dtype=np.int16);dc=dr.copy()
    for code,(dy,dx) in OFF.items():dr[code]=dy;dc[code]=dx
    nr=rr+dr[fd];nc=cc+dc[fd];dest=nr*W+nc
    candidate=np.isin(fd,list(OFF))&(nr>=0)&(nr<H)&(nc>=0)&(nc<W)
    pos=np.searchsorted(gid,dest);safe=np.clip(pos,0,N-1)
    found=candidate&(pos<N)&(gid[safe]==dest)
    down=np.where(found,pos,-1).astype(np.int32)
    bad=np.zeros(N,dtype=bool);valid=np.where(found)[0];bad[valid]=A[down[valid]]<A[valid]-1e-6
    if bad.any():raise ValueError(f'{bad.sum()} enlaces D8 con área decreciente aguas abajo.')
    step=np.hypot(dr[fd]*abs(tr.e),dc[fd]*abs(tr.a))
    indeg=np.bincount(down[down>=0],minlength=N).astype(np.int16);work=indeg.copy()
    queue=deque(np.flatnonzero(work==0).tolist());topo=np.empty(N,dtype=np.int32);j=0
    maxorder=np.ones(N,dtype=np.int16);countmax=np.zeros(N,dtype=np.int16)
    order=np.ones(N,dtype=np.int16);main=np.full(N,-1,dtype=np.int32)
    print('Nodos de red global:',N,'· ordenando conectividad',flush=True)
    while queue:
        u=queue.popleft();topo[j]=u;j+=1
        order[u]=maxorder[u]+int(countmax[u]>=2)
        v=down[u]
        if v>=0:
            if order[u]>maxorder[v]:maxorder[v]=order[u];countmax[v]=1
            elif order[u]==maxorder[v]:countmax[v]+=1
            if main[v]<0 or A[u]>A[main[v]] or (A[u]==A[main[v]] and gid[u]<gid[main[v]]):main[v]=u
            work[v]-=1
            if work[v]==0:queue.append(int(v))
    if j!=N:raise ValueError('La red D8 contiene un ciclo; no se calculan χ ni ksn.')
    distance=np.zeros(N);chi=np.zeros(N)
    for u in topo[::-1]:
        v=down[u]
        if v>=0:
            distance[u]=distance[v]+step[u]
            chi[u]=chi[v]+step[u]*.5*((A0/A[u])**theta+(A0/A[v])**theta)
    g={'gid':gid,'A':A,'z':Z,'fd':fd,'down':down,'main':main,'order':order,'indegree':indeg,
       'distance':distance,'chi':chi,'step':step,'width':np.array(W),'height':np.array(H),
       'transform':np.array(list(tr)[:6]),'theta':np.array(theta),'A0':np.array(A0)}
    np.savez(target,**g)
    info=signature|{'n_nodos':N,'n_fuentes':int((indeg==0).sum()),'n_salidas':int((down<0).sum()),
                    'convencion':'ESRI D8; continuidad aguas arriba por mayor área; Strahler sobre red A>=umbral',
                    'raster_shape':[H,W],'cache':str(target)}
    meta.write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf8')
    return g,info

def xy(g,idx):
    tr=g['transform'];r=g['gid'][idx]//int(g['width']);c=g['gid'][idx]%int(g['width'])
    return tr[2]+(c+.5)*tr[0],tr[5]+(r+.5)*tr[4]

def ventana_ols(s,chi,z,W=1000.,min_n=5,min_fraction=.45,theta=.4,A0=1e6):
    s=np.asarray(s,float);x=np.asarray(chi,float);z=np.asarray(z,float)
    if len(s)<2:return {k:np.full(len(s),np.nan) for k in ['ksn','R2','RMSE','n','span','parcial']}
    if not(np.diff(s)>0).all():raise ValueError('Distancia del perfil no creciente.')
    x=x-x[0];y=z-z[0]
    lo=np.searchsorted(s,s-W/2);hi=np.searchsorted(s,s+W/2,side='right');n=hi-lo
    sums=[]
    for v in [x,y,x*x,y*y,x*y]:
        cs=np.r_[0,np.cumsum(v)];sums.append(cs[hi]-cs[lo])
    sx,sy,sxx,syy,sxy=sums
    den=sxx-sx*sx/n;vary=syy-sy*sy/n;cov=sxy-sx*sy/n
    span=s[np.maximum(hi-1,0)]-s[lo]
    good=(n>=min_n)&(den>1e-10)&(span>=W*min_fraction)
    slope=np.divide(cov,den,out=np.full(len(n),np.nan),where=good);good &= slope>0
    ksn=np.where(good,slope*A0**theta,np.nan)
    r2=np.divide(cov**2,den*vary,out=np.full(len(n),np.nan),where=good&(vary>0))
    rmse=np.sqrt(np.maximum(0,np.divide(vary-np.divide(cov**2,den,out=np.zeros(len(n)),where=den>0),n)))
    return {'ksn':ksn,'R2':np.clip(r2,0,1),'RMSE':np.where(good,rmse,np.nan),'n':n,'span':span,'parcial':span<W-25.}

def calcular_ksn(g,poly,window=1000.,min_n=5,min_fraction=.45):
    N=len(g['gid']);allidx=np.arange(N);x,y=xy(g,allidx)
    inside=shapely.intersects_xy(poly,x,y)
    starts=np.flatnonzero(g['main']<0);fields={k:np.full(N,np.nan) for k in ['ksn','R2','RMSE','n','span','parcial']}
    npaths=0
    for h in starts:
        path=[];u=int(h)
        while u>=0:
            path.append(u);v=int(g['down'][u])
            if v<0 or g['main'][v]!=u:break
            u=v
        owned=np.array(path,dtype=np.int32)
        if not inside[owned].any():continue
        # Extensión aguas abajo a través de la confluencia para no cortar la ventana.
        path_full=path.copy();v=int(g['down'][path[-1]]);last=path[-1]
        while v>=0 and g['distance'][last]-g['distance'][v]<=window/2+25:
            path_full.append(v);v=int(g['down'][v])
        idx=np.array(path_full[::-1],dtype=np.int32)
        stats=ventana_ols(g['distance'][idx],g['chi'][idx],g['z'][idx],window,min_n,min_fraction,float(g['theta']),float(g['A0']))
        # Los nodos ajenos sólo aportan contexto; no se sobrescriben sus propias ventanas.
        for key in fields:fields[key][owned]=stats[key][::-1][:len(owned)]
        npaths+=1
    print('Red dentro del polígono:',int(inside.sum()),'nodos · ksn válidos:',int(np.isfinite(fields['ksn'][inside]).sum()),flush=True)
    return fields,inside,{'n_cauces_continuos_procesados':npaths,'ventana_m':window,'min_nodos':min_n,
                          'fraccion_minima':min_fraction,'metodo':'OLS z–χ con intercepto libre; ksn=Mχ A0^theta'}

def segmentos_aoi(g,poly,inside,fields):
    u=np.flatnonzero((g['down']>=0)&(inside|np.where(g['down']>=0,inside[np.maximum(g['down'],0)],False)))
    v=g['down'][u];x1,y1=xy(g,u);x2,y2=xy(g,v)
    geom=shapely.linestrings(np.stack([np.c_[x1,y1],np.c_[x2,y2]],axis=1))
    clipped=shapely.intersection(geom,poly);keep=~shapely.is_empty(clipped)
    u=u[keep];v=v[keep];clipped=clipped[keep]
    d={'node_u':u,'node_v':v,'gid_u':g['gid'][u],'area_km2':g['A'][u]/1e6,'orden':g['order'][u],
       'chi_m':g['chi'][u],**{k:a[u] for k,a in fields.items()}}
    return gpd.GeoDataFrame(d,geometry=clipped,crs=CRS)

def leer_fallas(path,poly,nombres=None):
    f=gpd.read_file(path)
    if f.crs is None:raise ValueError('Fallas sin CRS.')
    f=f.to_crs(CRS);f=f[f.intersects(poly)].copy()
    f['nombre_base']=f.falla.astype(str).str.replace('*','',regex=False).str.strip()
    if nombres is not None:f=f[f.nombre_base.isin(nombres)]
    return f.reset_index(drop=True)

def cruzar_fallas(edges,fallas):
    j=gpd.sjoin(edges,fallas[['f_id','falla','nombre_base','geometry']],how='inner',predicate='intersects')
    records=[]
    for _,row in j.iterrows():
        feature=fallas.loc[int(row.index_right)]
        inter=row.geometry.intersection(feature.geometry)
        if inter.geom_type=='Point':pts=[inter]
        elif inter.geom_type in ['MultiPoint','GeometryCollection']:pts=[p for p in inter.geoms if p.geom_type=='Point']
        else:continue # Una superposición lineal no constituye un cruce puntual inequívoco.
        for pt in pts:records.append({'node_u':int(row.node_u),'node_v':int(row.node_v),'f_id':int(feature.f_id),
            'falla':feature.falla,'nombre_base':feature.nombre_base,'orden_cruce':int(row.orden),
            'x_cruce':pt.x,'y_cruce':pt.y,'geometry':pt})
    if not records:return gpd.GeoDataFrame(columns=['node_u','node_v','f_id','falla','nombre_base','orden_cruce','x_cruce','y_cruce','geometry'],crs=CRS)
    r=gpd.GeoDataFrame(records,crs=CRS)
    r['qx']=(r.x_cruce*1000).round().astype('int64');r['qy']=(r.y_cruce*1000).round().astype('int64')
    return r.drop_duplicates(['f_id','qx','qy']).drop(columns=['qx','qy']).reset_index(drop=True)

def seleccionar_perfiles(g,crossings,min_order=2,min_length=2000.):
    profiles={};selection=[];group={};ids={}
    for i,c in crossings.iterrows():
        if c.orden_cruce<min_order:
            selection.append({'cruce_fila':i,'estado':'excluido por orden de cabecera','orden':c.orden_cruce});continue
        upstream=[];u=int(c.node_u)
        while u>=0:upstream.append(u);u=int(g['main'][u])
        head=upstream[-1]
        if head not in group:
            p=[];u=head
            while u>=0:p.append(u);u=int(g['down'][u])
            idx=np.array(p[::-1],dtype=np.int32)
            length=float(g['distance'][idx[-1]]-g['distance'][idx[0]])
            if length<min_length:
                selection.append({'cruce_fila':i,'estado':'excluido por longitud','longitud_m':length});continue
            pid=f'F{len(group)+1:03d}';group[head]=idx;ids[head]=pid
            x,y=xy(g,idx)
            profiles[pid]=pd.DataFrame({'node':idx,'gid':g['gid'][idx],'x_utm':x,'y_utm':y,
                'distancia_m':g['distance'][idx]-g['distance'][idx[0]],'chi_m':g['chi'][idx],
                'elevacion_m':g['z'][idx],'area_m2':g['A'][idx],'orden':g['order'][idx]})
        selection.append({'cruce_fila':i,'estado':'seleccionado','perfil_id':ids[head],
                          'orden':c.orden_cruce,'longitud_m':float(profiles[ids[head]].distancia_m.iloc[-1])})
    audit=pd.DataFrame(selection)
    if not audit.empty:
        crossings=crossings.copy();crossings['cruce_fila']=np.arange(len(crossings));crossings=crossings.merge(audit,on='cruce_fila',how='left')
    for pid,p in profiles.items():
        q=crossings[crossings.get('perfil_id',pd.Series(index=crossings.index,dtype=str)).eq(pid)]
        for i,c in q.iterrows():
            u,v=int(c.node_u),int(c.node_v);xu,yu=xy(g,np.array([u]));xv,yv=xy(g,np.array([v]))
            edge=LineString([(xu[0],yu[0]),(xv[0],yv[0])]);frac=edge.project(Point(c.x_cruce,c.y_cruce))/edge.length
            crossings.loc[i,'chi_cruce_m']=g['chi'][u]+frac*(g['chi'][v]-g['chi'][u])
            crossings.loc[i,'dist_cruce_m']=g['distance'][u]+frac*(g['distance'][v]-g['distance'][u])-g['distance'][int(p.node.iloc[0])]
            crossings.loc[i,'z_cruce_m']=g['z'][u]+frac*(g['z'][v]-g['z'][u])
    print('Cruces totales:',len(crossings),'· perfiles completos seleccionados:',len(profiles),flush=True)
    return profiles,crossings,audit

def asignar_litologia(profiles,path):
    geology=gpd.read_file(path).to_crs(32718)
    if not geology.is_valid.all():geology.geometry=geology.geometry.make_valid()
    rows=[]
    for pid,p in profiles.items():
        q=p.copy();q['perfil_id']=pid;q['fila_perfil']=np.arange(len(q));rows.append(q)
    if not rows:return profiles,geology
    allp=pd.concat(rows,ignore_index=True)
    # Un gid es un centro de celda nativo; sus atributos espaciales no dependen del perfil.
    pts=allp[['gid','x_utm','y_utm']].drop_duplicates('gid').reset_index(drop=True)
    pts=gpd.GeoDataFrame(pts,geometry=gpd.points_from_xy(pts.x_utm,pts.y_utm),crs=32718)
    cols=[c for c in ['codigo','grupo_lito','unidad','fuente','estado','uso_metr','color_hex','geometry'] if c in geology]
    joined=gpd.sjoin(pts,geology[cols],how='left',predicate='intersects')
    a=pd.DataFrame(index=pts.index)
    a['gid']=pts.gid;a['codigo_lito']='Sin cobertura';a['descripcion_lito']='Sin cobertura cartográfica'
    a['color_lito']='#8b949e';a['estado_lito']='sin cobertura';a['fuentes_lito']=''
    hit=joined[joined.index_right.notna()]
    if len(hit):
        count=hit.groupby(level=0,sort=False).codigo.nunique(dropna=True)
        first=hit[~hit.index.duplicated(keep='first')]
        src=hit.groupby(level=0,sort=False).fuente.agg(lambda v:' / '.join(sorted(v.dropna().astype(str).unique())))
        a.loc[src.index,'fuentes_lito']=src
        ids=count[count.eq(1)].index
        code=hit[hit.codigo.notna()];code=code[~code.index.duplicated(keep='first')]
        a.loc[ids,'codigo_lito']=code.loc[ids,'codigo'].astype(str)
        a.loc[ids,'descripcion_lito']=first.loc[ids,'grupo_lito'].astype(str)+' · '+first.loc[ids,'unidad'].astype(str)
        a.loc[ids,'color_lito']=first.loc[ids,'color_hex'].astype(str)
        a.loc[ids,'estado_lito']='intersección con cartografía, no validación de roca del lecho'
        ids=count[count.gt(1)].index
        if len(ids):
            labels=hit.loc[hit.index.isin(ids)].groupby(level=0,sort=False).codigo.agg(lambda v:' / '.join(sorted(v.dropna().astype(str).unique())))
            a.loc[ids,'codigo_lito']='Superposición';a.loc[ids,'descripcion_lito']='Unidades superpuestas: '+labels
            a.loc[ids,'color_lito']='#4d4d4d';a.loc[ids,'estado_lito']='ambiguo; no se selecciona una unidad'
    result=allp.merge(a,on='gid',how='left',validate='many_to_one',sort=False)
    out={pid:q.drop(columns=['perfil_id','fila_perfil']).reset_index(drop=True) for pid,q in result.groupby('perfil_id',sort=False)}
    return out,geology


def cruces_en_perfiles(profiles,fallas):
    rows=[]
    for pid,p in profiles.items():
        line=LineString(p[['x_utm','y_utm']].to_numpy())
        s=np.r_[0,np.cumsum(np.hypot(np.diff(p.x_utm),np.diff(p.y_utm)))]
        for _,f in fallas[fallas.intersects(line)].iterrows():
            hit=line.intersection(f.geometry)
            points=[hit] if hit.geom_type=='Point' else [x for x in getattr(hit,'geoms',[]) if x.geom_type=='Point']
            for point in points:
                at=line.project(point)
                rows.append({'perfil_id':pid,'f_id':int(f.f_id),'falla':f.falla,'nombre_base':f.nombre_base,
                             'chi_cruce_m':float(np.interp(at,s,p.chi_m)),
                             'dist_cruce_m':float(np.interp(at,s,p.distancia_m)),
                             'z_cruce_m':float(np.interp(at,s,p.elevacion_m)),
                             'x_cruce':point.x,'y_cruce':point.y})
    return pd.DataFrame(rows)
