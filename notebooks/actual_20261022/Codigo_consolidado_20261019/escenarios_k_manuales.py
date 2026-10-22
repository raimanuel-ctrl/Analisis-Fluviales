"""Escenarios exhaustivos B × fuente de altura × cruce P_T; sin selección final."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import map_coordinates
from shapely.geometry import Point,LineString
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import revision_b_base as rb

CRS='EPSG:32718'
OFFSETS={1:(0,1),2:(1,1),4:(1,0),8:(1,-1),16:(0,-1),32:(-1,-1),64:(-1,0),128:(-1,1)}
COLORS={'5e':'#b30ba6','5c':'#138c45'}

def leer_puntos(path):
    d=pd.read_csv(path,keep_default_na=False,float_precision='round_trip')
    keys=['rid','caso_B','id_punto','etiqueta','x_utm18s_m','y_utm18s_m','elevacion_dem_m','crs']
    if not set(keys).issubset(d):raise ValueError('Faltan columnas del CSV de selección.')
    for col in ['rid','id_punto','x_utm18s_m','y_utm18s_m','elevacion_dem_m']:d[col]=pd.to_numeric(d[col],errors='raise')
    if not np.isfinite(d[['rid','id_punto','x_utm18s_m','y_utm18s_m','elevacion_dem_m']].to_numpy(float)).all():raise ValueError('Coordenadas o IDs no finitos.')
    if not d.rid.eq(d.rid.astype(int)).all() or not d.rid.between(1,18).all():raise ValueError('ID de río inválido.')
    if not d.crs.eq(CRS).all():raise ValueError('Se necesita EPSG:32718.')
    if d.duplicated(['caso_B','id_punto']).any():raise ValueError('ID de punto repetido dentro del mismo escenario de origen.')
    d['fila_CSV']=np.arange(2,len(d)+2)
    d['mis_contexto_visor']=d['mis']
    # Decisión explícita del usuario: la etiqueta prevalece sobre el contexto del visor.
    d['mis']=np.where(d.etiqueta.astype(str).str.strip().eq(''),'5e','5c')
    d['clasificacion_MIS']='sin etiqueta=5e; con etiqueta=5c, instrucción del usuario'
    d['id_B_corto']=[f'N{i:02d}' for i in range(1,len(d)+1)]
    d['caso_visor']=d.caso_B
    return d

def extender_perfil(base,p,rid,max_m=8000.):
    """Recupera el cauce aguas abajo con el D8 y área ya usados; no extrapola rectas."""
    base=Path(base);cfg=json.loads((base/'Configuracion.json').read_text(encoding='utf-8-sig'))
    topo=base/f"Productos/Datos/{cfg['run_id']}/topografia"
    files={'flow':topo/f"flow_direction_D8_12m_{cfg['run_id']}.tif",
           'area':topo/f"area_aportante_12m_{cfg['run_id']}.tif",'dem':base/cfg['dem']}
    mf=json.loads((topo/f"manifest_{cfg['run_id']}.json").read_text(encoding='utf8'))
    for key,name in [('flow','flow_direction'),('area','area')]:
        if rb.huella(files[key])!=mf['files'][name]['sha256']:raise ValueError(f'Cambió el raster {key}.')
    x0,y0=p[['x_utm','y_utm']].iloc[0];arrays={};transform=None
    for key,path in files.items():
        with rasterio.open(path) as src:
            if src.crs!=rasterio.crs.CRS.from_epsg(32718):raise ValueError('CRS de raster distinto de UTM18S.')
            if transform is None:
                transform=src.transform
                win=from_bounds(x0-max_m-100,y0-max_m-100,x0+max_m+100,y0+max_m+100,transform).round_offsets().round_lengths()
                tr=src.window_transform(win)
            elif src.transform!=transform:raise ValueError('DEM, D8 y área no están alineados.')
            arrays[key]=src.read(1,window=win,masked=True,boundless=True).astype(float).filled(np.nan)
    inv=~tr
    def cell(x,y):
        c,r=inv*(x,y);return int(np.floor(r)),int(np.floor(c))
    i,j=cell(x0,y0)
    if not np.isclose(arrays['area'][i,j],p.area_m2.iloc[0],atol=1e-5):raise ValueError('Área del outlet no coincide con el perfil vigente.')
    if not np.isclose(arrays['dem'][i,j],p.elevacion_m.iloc[0],atol=1e-5):raise ValueError('Z del outlet no coincide con el perfil vigente.')
    # Comprobar convención de flujo con nodos consecutivos ya conocidos del río.
    for k in range(1,min(60,len(p))):
        a,b=cell(p.x_utm.iloc[k],p.y_utm.iloc[k]);di,dj=OFFSETS.get(int(arrays['flow'][a,b]),(999,999))
        if (a+di,b+dj)!=cell(p.x_utm.iloc[k-1],p.y_utm.iloc[k-1]):raise ValueError('D8 no reproduce el perfil vigente.')
    trail=[];seen=set();length=0.;stop=''
    while True:
        if (i,j) in seen:raise ValueError('Bucle D8 aguas abajo.')
        if not(0<=i<arrays['flow'].shape[0] and 0<=j<arrays['flow'].shape[1]):stop='límite de reserva';break
        seen.add((i,j));x,y=tr*(j+.5,i+.5);z=arrays['dem'][i,j];area=arrays['area'][i,j]
        if not np.isfinite([z,area]).all() or area<=0:stop='fin de datos válidos';break
        trail.append({'x_utm':x,'y_utm':y,'elevacion_m':z,'area_m2':area})
        value=arrays['flow'][i,j]
        if not np.isfinite(value) or value not in OFFSETS:stop='fin de flujo D8';break
        di,dj=OFFSETS[int(value)];step=np.hypot(di*tr.e,dj*tr.a)
        if length+step>max_m:stop='límite de longitud';break
        length+=step;i+=di;j+=dj
    if len(trail)<2:raise ValueError('No se recuperó cauce aguas abajo.')
    chain=pd.DataFrame(trail)
    if (np.diff(chain.area_m2)<-1e-6).any():raise ValueError('El área disminuye aguas abajo en la extensión.')
    ext=chain.iloc[:0:-1].copy()
    prefix_len=float(np.hypot(np.diff(chain.x_utm),np.diff(chain.y_utm)).sum())
    ext['distancia_m']=np.r_[0,np.cumsum(np.hypot(np.diff(ext.x_utm),np.diff(ext.y_utm)))]
    old=p[['x_utm','y_utm','elevacion_m','area_m2','distancia_m']].copy();old.distancia_m+=prefix_len
    ext['tramo_recuperado_aguas_abajo']=True;old['tramo_recuperado_aguas_abajo']=False
    full=rb.preparar_perfil(pd.concat([ext,old],ignore_index=True))
    shift=full.chi_m.iloc[len(ext)]
    error=float(np.max(np.abs(full.chi_m.iloc[len(ext):].to_numpy()-shift-p.chi_m.to_numpy())))
    if error>1e-7:raise ValueError('La extensión alteró diferencias χ del perfil original.')
    return full,{'rid':rid,'n_nodos_recuperados':len(ext),'longitud_recuperada_m':prefix_len,
        'desplazamiento_chi_m':shift,'max_error_delta_chi_original_m':error,'fin_trazado':stop,
        'nota':'D8 y área vigentes; el tramo recuperado puede ser aluvial. No se valida incisión en roca en él.',
        'fuentes':{str(file):rb.huella(file) for file in files.values()}}

def muestra_dem(path,xy):
    with rasterio.open(path) as src:
        x,y=xy;win=from_bounds(x-36,y-36,x+36,y+36,src.transform).round_offsets().round_lengths()
        a=src.read(1,window=win,masked=True,boundless=True).astype(float).filled(np.nan);inv=~src.window_transform(win)
        c,r=inv*(x,y)
        return float(map_coordinates(a,[[r-.5],[c-.5]],order=1,mode='constant',cval=np.nan,prefilter=False)[0])

def construir_casos(base,perfiles,puntos,ruta_marcas,exclusiones,dem,max_separacion_m=100.):
    marks=gpd.read_file(ruta_marcas)
    excluded=set(json.loads(Path(exclusiones).read_text(encoding='utf8'))['sha256_wkb'])
    fids=[int(fid) for fid,row in marks.iterrows() if hashlib.sha256(row.geometry.wkb).hexdigest() in excluded]
    old,qc=rb.construir_B(perfiles,marks,None,excluir_fid=fids)
    old=old.sort_values(['rid','mis','caso_B']).reset_index(drop=True)
    old['origen_B']='anterior';old['id_B_corto']=[f'A{i:02d}' for i in range(1,len(old)+1)]
    old['x_seleccion']=old.x_B;old['y_seleccion']=old.y_B
    old['separacion_al_cauce_m']=0.;old['proyectado_en_extremo']=False
    old['H_dem_punto_m']=[muestra_dem(dem,(x,y)) for x,y in zip(old.x_B,old.y_B)]
    old['estado_B']='interpretación anterior; no aceptada automáticamente'
    old['B_computable']=True
    new=[]
    for _,v in puntos.iterrows():
        rid=int(v.rid);p=perfiles[rid];line=rb.linea_rio(p);point=Point(v.x_utm18s_m,v.y_utm18s_m)
        along=line.project(point);at=line.interpolate(along);dist=point.distance(at)
        endpoint=along<1e-6 or along>line.length-1e-6
        usable=not endpoint and dist<=max_separacion_m
        values={k:float(np.interp(along,p.s_geom_m,p[col])) for k,col in [('s_B_m','distancia_m'),('chi_B_m','chi_m'),('z_B_m','elevacion_m')]}
        new.append({'caso_B':f'R{rid:02d}_{v.mis}_nuevo_{v.id_B_corto}','rid':rid,'mis':v.mis,
          'grupo':'manual_3D','escenario':v.id_B_corto,'origen_B':'nuevo 3D','id_B_corto':v.id_B_corto,
          'tipo_B':'quiebre seleccionado en visor, proyectado al cauce','revisado':False,
          'x_seleccion':v.x_utm18s_m,'y_seleccion':v.y_utm18s_m,'x_B':at.x,'y_B':at.y,
          'H_dem_punto_m':v.elevacion_dem_m,'separacion_al_cauce_m':dist,'proyectado_en_extremo':endpoint,
          'B_computable':usable,'estado_B':('proyección local; revisar geometría' if usable else 'sin B utilizable: fuera del cauce/perfil'),
          'fila_CSV':v.fila_CSV,'id_punto_visor':v.id_punto,'etiqueta_visor':v.etiqueta,'mis_contexto_visor':v.mis_contexto_visor,
          'caso_visor':v.caso_visor,'clasificacion_MIS':v.clasificacion_MIS,
          'fuente':str(Path(base)/'Fuentes/Puntos_B_todos_los_rios.csv'),**values})
    cases=pd.concat([old,pd.DataFrame(new)],ignore_index=True)
    cases['B_en_extension_recuperada']=[float(np.interp(b.s_B_m,perfiles[int(b.rid)].distancia_m,
        perfiles[int(b.rid)]['tramo_recuperado_aguas_abajo'].astype(float)))>.5 for _,b in cases.iterrows()]
    # Los B inválidos conservan XY de diagnóstico; jamás reciben χ calculable.
    cases.loc[~cases.B_computable,['chi_B_m','s_B_m','z_B_m']]=np.nan
    return cases,qc

def asociar_alturas(cases,terrazas,nombres):
    t=terrazas.copy();t['mis']=t.mis_asignado.map(rb.mis_normal)
    rows=[]
    for _,b in cases.iterrows():
        same=t[t.mis.eq(b.mis)].copy()
        if same.empty:raise ValueError('Sin referencia J-M del MIS solicitado.')
        same['dist_B']=np.hypot(same.este_utm18s_m-b.x_B,same.norte_utm18s_m-b.y_B)
        v=same.sort_values(['dist_B','terrace_id']).iloc[0]
        common=b.to_dict()|{'rio':nombres[int(b.rid)],'terraza_id':v.terrace_id,'H_JM_m':v.elevacion_m,
          'error_H_JM_m':v.error_elevacion_m,'T_ka':v.edad_referencia_ka,'error_T_ka':v.incertidumbre_edad_ka,
          'x_terraza':v.este_utm18s_m,'y_terraza':v.norte_utm18s_m,'distancia_terraza_B_m':v.dist_B,
          'criterio_terraza':'mínima distancia euclidiana UTM a B sobre el cauce, restringida al mismo MIS',
          'edad_usada':'edad asignada al MIS, no edad directa de sedimento','fuente_H_JM':v.fuente}
        for fuente,H,err in [('JM',v.elevacion_m,v.error_elevacion_m),('DEM_punto',b.H_dem_punto_m,0.)]:
            rows.append(common|{'fuente_H':fuente,'H_m':H,'error_H_m':err,
                'caso_altura':f'{b.caso_B}_H_{fuente}',
                'nota_H':('altura publicada J-M' if fuente=='JM' else 'altura del punto seleccionado, fija; variante experimental, no validada como altura de terraza')})
    return pd.DataFrame(rows)

def calcular(refs,perfiles):
    records=[];qc=[]
    for _,b in refs.iterrows():
        if not b.B_computable or not np.isfinite(b.H_m):
            qc.append({'caso_altura':b.caso_altura,'n_cruces':0,'estado':'sin B/H utilizable'});continue
        p=perfiles[int(b.rid)];cr=rb.cruces_cota(p,b.H_m)
        # También conservar una meseta que comienza antes de B y termina después.
        cr=cr[(cr.s_PT_m>b.s_B_m+1e-6)|(cr.direccion.eq('meseta') & (cr.s_PT_fin_m>b.s_B_m+1e-6))]
        qc.append({'caso_altura':b.caso_altura,'n_cruces':len(cr),'n_ascendentes':int(cr.direccion.eq('ascendente').sum()),
                   'H_m':b.H_m,'min_Z_perfil':p.elevacion_m.min(),'max_Z_perfil':p.elevacion_m.max(),
                   'estado':'con cruces' if len(cr) else 'sin cruce de H aguas arriba de B'})
        for _,pt in cr.iterrows():
            flat=pt.direccion=='meseta';delta=pt.chi_PT_m-b.chi_B_m;den=rb.A0**rb.THETA*b.T_ka*1000
            K=delta/den if not flat and delta>0 else np.nan
            start=max(pt.s_PT_m,b.s_B_m)
            lower=(np.interp(start,p.distancia_m,p.chi_m)-b.chi_B_m)/den if flat else K
            upper=(pt.chi_PT_fin_m-b.chi_B_m)/den if flat else K
            eH=(b.error_H_m/abs(pt.Mchi_local)/den) if not flat and pt.Mchi_local!=0 else np.nan
            eT=K*b.error_T_ka/b.T_ka if np.isfinite(K) else np.nan
            approx=np.hypot(eH,eT) if np.isfinite(eH) else np.nan
            records.append(b.to_dict()|pt.to_dict()|{'escenario_K':f'{b.caso_altura}_PT{int(pt.cruce_PT):03d}',
                'delta_chi_m':delta,'K_m02_anio':K,'K_intervalo_min':lower,'K_intervalo_max':upper,
                'error_K_H_aprox':eH,'error_K_T_aprox':eT,'error_K_propagado_aprox':approx,
                'nota_error':'magnitudes publicadas H y T; aproximación local, sin error geométrico ni de datum; DEM fijo',
                'seleccionado_final':False,'estado':'escenario condicional, sin elección final',
                'datum_vertical_confirmado':False,'ascendente':pt.direccion=='ascendente',
                'PT_en_extension_recuperada':np.interp(pt.s_PT_m,p.distancia_m,p.tramo_recuperado_aguas_abajo.astype(float))>.5})
    return pd.DataFrame(records),pd.DataFrame(qc)

def resumir(resultados,refs):
    out=[]
    for _,b in refs.iterrows():
        q=resultados[resultados.caso_altura.eq(b.caso_altura)] if not resultados.empty else resultados
        for selection in ['todos_puntuales','solo_ascendentes']:
            x=q[q.K_m02_anio.gt(0)] if not q.empty else q
            if selection=='solo_ascendentes' and not x.empty:x=x[x.ascendente]
            k=x.K_m02_anio if not x.empty else pd.Series(dtype=float)
            med=k.median()
            out.append({'rid':b.rid,'rio':b.rio,'caso_B':b.caso_B,'id_B':b.id_B_corto,'mis':b.mis,'origen_B':b.origen_B,
                'fuente_H':b.fuente_H,'H_m':b.H_m,'terraza_id':b.terraza_id,'seleccion_cruces':selection,
                'n_K':len(k),'K_min':k.min(),'K_mediana':med,'K_max':k.max(),
                'rango_relativo_pct':100*(k.max()-k.min())/med if len(k) and med>0 else np.nan,
                'nota':'mínimo/mediana/máximo descriptivos de escenarios, no estimador aceptado ni intervalo de confianza'})
    return pd.DataFrame(out)

def graficos_por_rio(resultados,refs,perfiles,dem,folder,escala='linear',mostrar_descendentes=True):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    files=[]
    for rid,g in refs.groupby('rid'):
        r=resultados[resultados.rid.eq(rid)].copy() if not resultados.empty else resultados
        if r.empty:continue
        order=g.drop_duplicates('caso_B').sort_values(['origen_B','mis','caso_B'])
        ids=order.caso_B.tolist();pos={v:i for i,v in enumerate(ids)}
        fig,axs=plt.subplots(2,1,figsize=(max(10,len(ids)*1.5),8.4),sharex=True,sharey=True,layout='constrained')
        for ax,source in zip(axs,['JM','DEM_punto']):
            h=r[r.fuente_H.eq(source)]
            for case,z in h.groupby('caso_B'):
                k=z[z.K_m02_anio.gt(0)];i=pos[case];color=COLORS[z.mis.iloc[0]]
                if not k.empty:ax.plot([i,i],[k.K_m02_anio.min(),k.K_m02_anio.max()],color='#737c86',lw=2,zorder=1)
                for direction,marker in [('ascendente','o'),('descendente','v')]:
                    if direction=='descendente' and not mostrar_descendentes:continue
                    q=k[k.direccion.eq(direction)].sort_values('s_PT_m')
                    jitter=np.linspace(-.14,.14,len(q)) if len(q)>1 else np.zeros(len(q))
                    ax.scatter(i+jitter,q.K_m02_anio,s=35,marker=marker,
                        facecolors=color if direction=='ascendente' else 'none',edgecolors=color,zorder=3)
                for _,flat in z[z.direccion.eq('meseta')].iterrows():
                    ax.plot([i+.23,i+.23],[flat.K_intervalo_min,flat.K_intervalo_max],color=color,lw=4,alpha=.4)
            ax.set_title('Altura publicada de la terraza J-M más cercana' if source=='JM' else 'Altura DEM del punto marcado · experimento',loc='left',fontsize=12)
            ax.set_ylabel(r'$K$ (m$^{0.2}$ año$^{-1}$)');ax.grid(axis='y',alpha=.2);ax.set_yscale(escala)
            if escala=='linear':
                formatter=ScalarFormatter(useMathText=True);formatter.set_powerlimits((0,0));ax.yaxis.set_major_formatter(formatter);ax.set_ylim(bottom=0)
        axs[-1].set_xticks(range(len(ids)),[f'{row.id_B_corto} · {row.mis}\n'+('B anterior' if row.origen_B=='anterior' else 'B nuevo 3D') for _,row in order.iterrows()])
        from matplotlib.lines import Line2D
        legend=[Line2D([],[],marker='o',ls='',color=COLORS[m],label=f'MIS {m}') for m in ['5e','5c']]
        legend += [Line2D([],[],marker='o',ls='',color='#333',label='P_T ascendente'),Line2D([],[],marker='v',ls='',mfc='none',color='#333',label='P_T descendente'),Line2D([],[],color='#737c86',label='Rango de escenarios, no error')]
        fig.legend(handles=legend,loc='outside lower center',ncol=3,fontsize=9)
        fig.suptitle(f'R{int(rid):02d} · {g.rio.iloc[0]} — variación de K por B y P_T\nCada símbolo es un cruce; ningún escenario ha sido seleccionado',fontsize=14)
        path=folder/f'R{int(rid):02d}_K_por_B_y_altura.png';fig.savefig(path,dpi=170,facecolor='white');plt.close(fig);files.append(path)
        # Perfil χ: puntos de partida B, cota de cada variante y todos los P_T.
        p=perfiles[int(rid)];fig,axs=plt.subplots(2,1,figsize=(12,8),layout='constrained')
        for ax,source in zip(axs,['JM','DEM_punto']):
            hh=r[r.fuente_H.eq(source)]
            maxchi=max(float(hh.chi_PT_m.max()) if len(hh) else 0,float(g.chi_B_m.max()))
            sub=p[p.chi_m<=maxchi+max(100,.1*maxchi)]
            ax.plot(sub.chi_m,sub.elevacion_m,color='#344b5c',lw=1.2,label='Perfil del cauce')
            for _,b in g[g.fuente_H.eq(source)].iterrows():
                color=COLORS[b.mis];ax.axhline(b.H_m,color=color,alpha=.35,lw=.8,ls='--')
                ax.plot([b.chi_B_m],[b.z_B_m],marker='s' if b.origen_B=='anterior' else '*',color=color,ms=7)
                ax.annotate(b.id_B_corto,(b.chi_B_m,b.z_B_m),xytext=(4,6),textcoords='offset points',fontsize=8,color=color)
            for direction,marker in [('ascendente','o'),('descendente','v')]:
                q=hh[hh.direccion.eq(direction)];ax.scatter(q.chi_PT_m,q.H_m,s=25,marker=marker,c=[COLORS[m] for m in q.mis],alpha=.6)
            ax.set_title(source+' · líneas discontinuas = alturas; círculos/triángulos = cruces P_T',loc='left',fontsize=11)
            ax.set_xlabel(r'$\chi$ (m), origen común del perfil de esta corrida');ax.set_ylabel('Elevación (m)');ax.grid(alpha=.15)
        fig.suptitle(f'R{int(rid):02d} · {g.rio.iloc[0]} — procedencia de Δχ\nΔχ = χ(P_T) − χ(B); no se usa la diferencia de elevaciones para calcular K',fontsize=13)
        path=folder/f'R{int(rid):02d}_perfil_chi_B_PT.png';fig.savefig(path,dpi=170,facecolor='white');plt.close(fig);files.append(path)
    return files

def tablas_png(resultados,refs,folder,max_rows=24):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);files=[]
    for rid,g in refs.groupby('rid'):
        r=resultados[resultados.rid.eq(rid)].copy() if not resultados.empty else resultados
        if r.empty:continue
        r=r.sort_values(['caso_B','fuente_H','s_PT_m'])
        for part,start in enumerate(range(0,len(r),max_rows),1):
            q=r.iloc[start:start+max_rows];cells=[]
            for _,v in q.iterrows():
                kval=f'{v.K_m02_anio:.3e}' if np.isfinite(v.K_m02_anio) else f'[{v.K_intervalo_min:.2e}, {v.K_intervalo_max:.2e}]'
                cells.append([v.id_B_corto,v.mis,'J-M' if v.fuente_H=='JM' else 'DEM',f'{v.H_m:.2f}',v.terraza_id,
                    f'PT{int(v.cruce_PT):03d}',v.direccion,f'{v.chi_B_m:.2f}',f'{v.chi_PT_m:.2f}',f'{v.delta_chi_m:.2f}',kval])
            fig,ax=plt.subplots(figsize=(17,max(2.5,.34*len(cells)+1.7)));ax.axis('off')
            table=ax.table(cellText=cells,colLabels=['B','MIS','Altura\nfuente','H (m)','Terraza J-M','Cruce','Dirección', 'χ(B) (m)','χ(P_T) (m)','Δχ (m)','K (m^0.2/año)'],
                cellLoc='center',loc='center',colWidths=[.04,.04,.055,.06,.1,.06,.09,.08,.08,.08,.18])
            table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,1.6)
            for (i,j),cell in table.get_celld().items():
                cell.set_facecolor('#dce8f2' if i==0 else ('white' if i%2 else '#f0f4f7'))
                cell.set_text_props(color='#172d3e',weight='bold' if i==0 else 'normal');cell.set_edgecolor('#cbd5df')
            ax.set_title(f'R{int(rid):02d} · {g.rio.iloc[0]} — todos los escenarios, parte {part}\nA = B anterior; N = B nuevo. Altura DEM experimental. Sin elección final.',fontsize=13,pad=15)
            path=folder/f'R{int(rid):02d}_tabla_K_{part:02d}.png';fig.savefig(path,dpi=160,bbox_inches='tight',facecolor='white');plt.close(fig);files.append(path)
    return files
