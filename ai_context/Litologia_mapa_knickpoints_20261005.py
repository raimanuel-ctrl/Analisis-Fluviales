# Code-only reading snapshot; not an independently executable script.
# Notebook is authoritative. Cell magics, external data and cell order still apply.

# Notebook cell 2
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json,hashlib
import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.integrate import cumulative_trapezoid
from shapely.geometry import Point
from IPython.display import display

BASE=Path(r'C:\datos\Memoria\Codigos Stanford')
RUN='20261005_post_auditoria_01'
PRODUCTOS=BASE/'Productos'/'post_auditoria'/RUN
GEO_DIR=Path(r'C:\datos\Memoria_Codex\outputs\Geologia_integrada_Chanco_Jorquera_Gajardo')
CRS='EPSG:32718'
RIOS={'3':'Estero Quebrada Honda','4':'Estero Malhueco','7':'Estero de Chanco (prueba)',
      '8':'Río Rahue','9':'Río Curanilahue','10':'Río Curanipe'}
RIOS_VER=list(RIOS)  # Para probar una figura: RIOS_VER=['9']
THETA,A0=0.40,1e6
INCLUIR_CANDIDATOS=True  # Candidatos del detector, NO knickpoints validados.
GUARDAR_TABLAS=True
GUARDAR_FIGURAS=False
FECHA=datetime.now(ZoneInfo('America/Los_Angeles')).strftime('%Y%m%d_%H%M%S_%f')
SALIDA=BASE/'Productos'/'litologia_y_mapa'/FECHA
entradas=[]
def registrar(path):
    path=Path(path)
    entradas.append(dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
def leer_manifest(etapa):
    path=PRODUCTOS/etapa/f'manifest_{RUN}.json'
    raw=path.read_bytes()
    entradas.append(dict(path=str(path),sha256=hashlib.sha256(raw).hexdigest()))
    m=json.loads(raw.decode('utf-8-sig'))
    if m.get('estado')!='completo': raise ValueError(f'Etapa incompleta: {etapa}')
    return m
def leer_tabla(m,key):
    item=m['files'][key]; p=Path(item['path'])
    registrar(p)
    if entradas[-1]['sha256']!=item['sha256']: raise ValueError(f'Archivo modificado: {p.name}')
    return pd.read_csv(p,dtype={'ID':str},float_precision='round_trip')

manifest=leer_manifest('perfiles')
if manifest['config'].get('crs')!=CRS:
    raise ValueError('Verificar CRS de los perfiles antes de cruzar con cartografía.')
perfiles={}
for rid in RIOS:
    p=leer_tabla(manifest,f'rio_{rid}')
    s,z,a=(p[c].to_numpy(float) for c in ['distancia_m','elevacion_m','area_m2'])
    if len(p)<2 or not np.isfinite(s).all() or not np.isclose(s[0],0) or np.any(np.diff(s)<=0):
        raise ValueError(f'{rid}: distancia inválida.')
    if not np.isfinite(a).all() or np.any(a<=0) or not np.isfinite(p[['x_utm','y_utm']]).all().all():
        raise ValueError(f'{rid}: área/coordenadas inválidas.')
    p['chi_040_m']=cumulative_trapezoid((A0/a)**THETA,x=s,initial=0)
    perfiles[rid]=p

for nombre in ['Geologia_unidades','Geologia_fallas']:
    for ext in ['.shp','.shx','.dbf','.prj','.cpg']:
        path=GEO_DIR/(nombre+ext)
        if path.exists(): registrar(path)
geo=gpd.read_file(GEO_DIR/'Geologia_unidades.shp')
fallas=gpd.read_file(GEO_DIR/'Geologia_fallas.shp')
if geo.crs is None or fallas.crs is None: raise ValueError('Cartografía sin CRS definido.')
geo=geo.to_crs(CRS); fallas=fallas.to_crs(CRS)
requeridas={'codigo','grupo_lito','fuente','color_hex','uso_metr','escala','unidad'}
if requeridas-set(geo): raise ValueError(f'Faltan campos: {requeridas-set(geo)}')
if geo.geometry.isna().any() or geo.geometry.is_empty.any() or not geo.geometry.is_valid.all():
    raise ValueError('Geología con geometrías nulas/vacías/inválidas; revisar fuente sin corregirla automáticamente.')
if fallas.geometry.isna().any() or fallas.geometry.is_empty.any() or not fallas.geometry.is_valid.all():
    raise ValueError('Fallas con geometrías inválidas; revisar fuente.')
display(geo[['fuente','escala','uso_metr']].drop_duplicates())
print('Perfiles cargados:',list(perfiles))


# Notebook cell 4
def asignar_unidades(puntos):
    # Un punto con dos polígonos sigue siendo ambiguo incluso si comparten código.
    joined=gpd.sjoin(puntos[['geometry']],geo[list(requeridas)+['geometry']],how='left',predicate='intersects')
    rows=[]
    for idx,g in joined.groupby(level=0,sort=False):
        valid=g[g.index_right.notna()]
        if len(valid)==1:
            r=valid.iloc[0]
            row={k:r[k] for k in requeridas}
            row.update(estado_geo='asignado',etiqueta=f"{r['codigo']} · {r['fuente']}")
        else:
            label='Sin cobertura' if len(valid)==0 else 'Ambiguo'
            row={k:None for k in requeridas}
            row.update(estado_geo=label,etiqueta=label,color_hex='#B7B7B7')
        rows.append(dict(indice=idx,**row))
    return pd.DataFrame(rows).set_index('indice').reindex(puntos.index)

segmentos={}; cambios=[]; propuestas=[]
for rid,p in perfiles.items():
    q=pd.DataFrame({'s0':p.distancia_m.to_numpy()[:-1],'s1':p.distancia_m.to_numpy()[1:],
                    'z0':p.elevacion_m.to_numpy()[:-1],'z1':p.elevacion_m.to_numpy()[1:],
                    'chi0':p.chi_040_m.to_numpy()[:-1],'chi1':p.chi_040_m.to_numpy()[1:]})
    puntos=gpd.GeoDataFrame(geometry=gpd.points_from_xy(
        (p.x_utm.to_numpy()[:-1]+p.x_utm.to_numpy()[1:])/2,
        (p.y_utm.to_numpy()[:-1]+p.y_utm.to_numpy()[1:])/2),crs=CRS)
    q=pd.concat([q,asignar_unidades(puntos)],axis=1)
    q['ID']=rid
    segmentos[rid]=q
    for j in range(1,len(q)):
        a,b=q.iloc[j-1],q.iloc[j]
        if a.etiqueta!=b.etiqueta:
            tipo=('cobertura_o_solape' if a.estado_geo!='asignado' or b.estado_geo!='asignado'
                  else 'cambio_fuente_posible_costura' if a.fuente!=b.fuente else 'cambio_unidad_muestreado')
            cambios.append(dict(ID=rid,distancia_m=float(b.s0),antes=a.etiqueta,despues=b.etiqueta,tipo=tipo))
    aluv=q.estado_geo.eq('asignado')&q.grupo_lito.fillna('').str.contains('fluvial|aluvial',case=False,regex=True)
    # Conservar unidades/fuentes separadas aunque dos intervalos se toquen.
    grupos=(aluv.ne(aluv.shift())|q.etiqueta.ne(q.etiqueta.shift())).cumsum()
    for _,g in q[aluv].groupby(grupos[aluv]):
        propuestas.append(dict(ID=rid,inicio_m=float(g.s0.iloc[0]),fin_m=float(g.s1.iloc[-1]),
            unidad=g.etiqueta.iloc[0],grupo=g.grupo_lito.iloc[0],estado='propuesta_no_aplicada'))
cambios_geo=pd.DataFrame(cambios,columns=['ID','distancia_m','antes','despues','tipo'])
propuestas_aluviales=pd.DataFrame(propuestas,columns=['ID','inicio_m','fin_m','unidad','grupo','estado'])
display(propuestas_aluviales)
display(pd.concat(segmentos.values()).groupby(['ID','estado_geo']).size().rename('segmentos').reset_index())


# Notebook cell 6
aceptados=[]; rechazados=[]; tablas=[]
ruta_analisis=PRODUCTOS/'analisis_perfiles'/f'manifest_{RUN}.json'
if ruta_analisis.exists():
    ma=leer_manifest('analisis_perfiles')
    tablas.append(('manual',leer_tabla(ma,'manuales')))
    if INCLUIR_CANDIDATOS:
        tablas.append(('candidato',leer_tabla(ma,'catalogo_candidatos')))
else:
    print('No existe manifiesto de análisis: se dibujarán cauces y geología sin puntos KP.')
for origen,tabla in tablas:
    tabla=tabla[tabla.ID.isin(RIOS)]
    for _,r in tabla.iterrows():
        rid=str(r.ID)
        try:
            p=perfiles[rid]; i=float(r['indice'])
            if not np.isfinite(i) or not i.is_integer() or not 0<=i<len(p): raise ValueError('Índice inválido')
            i=int(i); nodo=p.iloc[i]
            if r['perfil_hash']!=manifest['rios'][rid]['perfil_hash']: raise ValueError('Huella incompatible')
            for col in ['distancia_m','x_utm','y_utm','elevacion_m']:
                if not np.isclose(float(r[col]),float(nodo[col]),rtol=0,atol=.01): raise ValueError(f'Desajuste {col}')
            tipo='manual' if origen=='manual' or str(r.get('metodo','')).lower()=='manual' else 'candidato'
            aceptados.append(dict(ID=rid,id_kp=str(r.get('id_kp','')),tipo=tipo,metodo=r.get('metodo',''),
                theta_origen=r.get('theta',r.get('theta_seleccion',np.nan)),indice=i,
                distancia_m=float(nodo.distancia_m),elevacion_m=float(nodo.elevacion_m),
                chi_040_m=float(nodo.chi_040_m),x_utm=float(nodo.x_utm),y_utm=float(nodo.y_utm)))
        except (KeyError,ValueError,TypeError) as exc:
            rechazados.append(dict(ID=rid,id_kp=r.get('id_kp',''),origen=origen,motivo=str(exc)))
cols=['ID','id_kp','tipo','metodo','theta_origen','indice','distancia_m','elevacion_m','chi_040_m','x_utm','y_utm']
kp=pd.DataFrame(aceptados,columns=cols).drop_duplicates(['ID','id_kp','tipo']).reset_index(drop=True)
kp_rechazados=pd.DataFrame(rechazados,columns=['ID','id_kp','origen','motivo'])
if not kp.empty:
    puntos=gpd.GeoDataFrame(geometry=gpd.points_from_xy(kp.x_utm,kp.y_utm),crs=CRS)
    kp=pd.concat([kp,asignar_unidades(puntos)],axis=1)
print('Manuales compatibles:',int(kp.tipo.eq('manual').sum()),'; candidatos:',int(kp.tipo.eq('candidato').sum()))
display(kp)
display(kp_rechazados)


# Notebook cell 8
if GUARDAR_TABLAS:
    SALIDA.mkdir(parents=True,exist_ok=True)
    for nombre,tabla in [('segmentos_litologia',pd.concat(segmentos.values(),ignore_index=True)),
        ('cambios_cartograficos',cambios_geo),('propuestas_aluviales',propuestas_aluviales),
        ('kp_litologia',kp),('kp_rechazados',kp_rechazados)]:
        tabla.to_csv(SALIDA/f'{nombre}_{FECHA}.csv',index=False,encoding='utf-8-sig')
    (SALIDA/f'procedencia_{FECHA}.json').write_text(json.dumps(dict(
        fecha=FECHA,crs=CRS,rios=RIOS,theta=THETA,A0=A0,entradas=entradas,
        metodo_litologia='punto medio de segmentos; solapes ambiguos',
        origen_candidatos='exportación previa; no se repite detección',
        uso='exploratorio; sin distancias métricas a fallas/contactos'),indent=2,ensure_ascii=False),encoding='utf-8')
    print('Tablas guardadas:',SALIDA)


# Notebook cell 10
get_ipython().run_line_magic('matplotlib','inline')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection,PolyCollection
from matplotlib.lines import Line2D
from matplotlib.colors import is_color_like

GROSOR=2.4
BORDE_NEGRO=True
RELLENO=False
BANDA=False
VERTICALES=False
ETIQUETAS_KP=True
GRILLA=True
# BORDE_NEGRO=False
# GROSOR=3.5
# RELLENO=True
# BANDA=True
# VERTICALES=True  # Puede recargar perfiles con muchos cambios.
# ETIQUETAS_KP=False
# RIOS_VER=['9']

def colores(q):
    return [v if isinstance(v,str) and is_color_like(v) else '#B7B7B7' for v in q.color_hex]
def guardar_fig(fig,nombre):
    if GUARDAR_FIGURAS:
        SALIDA.mkdir(parents=True,exist_ok=True)
        fig.savefig(SALIDA/f'{nombre}_{FECHA}.png',dpi=250,bbox_inches='tight')
def dibujar_perfil(ax,rid,eje='distancia',relativo=False):
    q=segmentos[rid]; p=perfiles[rid]
    x0=q.s0.to_numpy()/1000 if eje=='distancia' else q.chi0.to_numpy()/1000
    x1=q.s1.to_numpy()/1000 if eje=='distancia' else q.chi1.to_numpy()/1000
    ref=float(p.elevacion_m.iloc[0]) if relativo else 0.
    z0=q.z0.to_numpy()-ref; z1=q.z1.to_numpy()-ref
    valid=np.isfinite(z0)&np.isfinite(z1)
    seg=np.stack([np.c_[x0,z0],np.c_[x1,z1]],axis=1)
    c=np.array(colores(q))
    if BORDE_NEGRO: ax.add_collection(LineCollection(seg[valid],colors='black',linewidths=GROSOR+1.4))
    ax.add_collection(LineCollection(seg[valid],colors=c[valid],linewidths=GROSOR))
    baseline=float(p.elevacion_m.iloc[0])-ref
    if RELLENO:
        polys=[[(a,baseline),(a,z),(b,w),(b,baseline)] for a,b,z,w in zip(x0[valid],x1[valid],z0[valid],z1[valid])]
        ax.add_collection(PolyCollection(polys,facecolors=c[valid],edgecolors='none',alpha=.22))
    if BANDA:
        polys=[[(a,0),(a,.045),(b,.045),(b,0)] for a,b in zip(x0,x1)]
        ax.add_collection(PolyCollection(polys,facecolors=c,edgecolors='none',transform=ax.get_xaxis_transform()))
    if VERTICALES:
        for s in cambios_geo.loc[(cambios_geo.ID==rid)&(cambios_geo.tipo=='cambio_unidad_muestreado'),'distancia_m']:
            v=s/1000 if eje=='distancia' else np.interp(s,p.distancia_m,p.chi_040_m)/1000
            ax.axvline(v,color='gray',ls=':',lw=.6,alpha=.5)
    for tipo,marker in [('manual','o'),('candidato','^')]:
        puntos=kp[(kp.ID==rid)&(kp.tipo==tipo)]
        x=puntos.distancia_m/1000 if eje=='distancia' else puntos.chi_040_m/1000
        ax.scatter(x,puntos.elevacion_m-ref,marker=marker,facecolors='white',edgecolors='black',s=38,zorder=5)
        if ETIQUETAS_KP:
            for xx,(_,r) in zip(x,puntos.iterrows()):
                ax.annotate(f"{tipo[0].upper()}{int(r.indice)}",(xx,r.elevacion_m-ref),xytext=(3,5),textcoords='offset points',fontsize=6)
    ax.autoscale_view()
    ax.set(xlim=(0,None),xlabel='Distancia desde outlet (km)' if eje=='distancia' else 'χ desde outlet (km)',
           ylabel='Elevación sobre outlet (m)' if relativo else 'Elevación (m)')
    ax.grid(GRILLA,alpha=.15)

for rid in RIOS_VER:
    fig,ax=plt.subplots(1,2,figsize=(13,4.5))
    dibujar_perfil(ax[0],rid)
    dibujar_perfil(ax[1],rid,'chi')
    ley=segmentos[rid].drop_duplicates('etiqueta')
    handles=[Line2D([0],[0],color=c,lw=4,label=label) for c,label in zip(colores(ley),ley.etiqueta)]
    handles += [Line2D([0],[0],color='black',marker=m,markerfacecolor='white',ls='',label=t)
                for t,m in [('Manual','o'),('Candidato exploratorio','^')]]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,0),ncol=3,fontsize=7)
    fig.suptitle(f'{rid} · {RIOS[rid]} · θ=0,40')
    fig.tight_layout(rect=(0,.12+len(handles)//3*.025,1,.95))
    guardar_fig(fig,f'litologia_rio_{rid}')
    plt.show(); plt.close(fig)


# Notebook cell 12
fig,ax=plt.subplots(figsize=(10,10))
xy=np.concatenate([perfiles[r][['x_utm','y_utm']].to_numpy() for r in RIOS_VER])
x0,y0=xy.min(axis=0)-2000; x1,y1=xy.max(axis=0)+2000
g=geo.cx[x0:x1,y0:y1]
g.plot(ax=ax,color=colores(g),alpha=.28,edgecolor='none')
f=fallas.cx[x0:x1,y0:y1]
if not f.empty: f.plot(ax=ax,color='#7B3444',linewidth=.65,alpha=.65)
for rid in RIOS_VER:
    p=perfiles[rid]
    ax.plot(p.x_utm,p.y_utm,lw=1.5,label=RIOS[rid])
    ax.annotate(rid,(p.x_utm.iloc[0],p.y_utm.iloc[0]),xytext=(4,4),textcoords='offset points',fontsize=9)
for tipo,marker in [('manual','o'),('candidato','^')]:
    k=kp[kp.ID.isin(RIOS_VER)&kp.tipo.eq(tipo)]
    if not k.empty: ax.scatter(k.x_utm,k.y_utm,marker=marker,facecolors='white',edgecolors='black',s=38,label=tipo,zorder=5)
ax.plot([],[],color='#7B3444',lw=.8,label='Trazas cartográficas preliminares')
ax.set(xlim=(x0,x1),ylim=(y0,y1),xlabel='Este UTM 18S (m)',ylabel='Norte UTM 18S (m)',title='Ríos y puntos compatibles · coordenadas reales')
ax.set_aspect('equal'); ax.ticklabel_format(style='plain',useOffset=False)
ax.legend(fontsize=8,loc='upper left',bbox_to_anchor=(1,1))
fig.tight_layout(); guardar_fig(fig,'mapa_rios_candidatos'); plt.show(); plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(13,4.5))
for rid in RIOS_VER:
    p=perfiles[rid]
    z=p.elevacion_m-float(p.elevacion_m.iloc[0])
    ax[0].plot(p.distancia_m/1000,z,lw=1.4,label=RIOS[rid])
    ax[1].plot(p.chi_040_m/1000,z,lw=1.4,label=RIOS[rid])
for eje,label in zip(ax,['Distancia desde outlet (km)','χ desde outlet (km)']):
    eje.set(xlim=(0,None),xlabel=label,ylabel='Elevación sobre outlet (m)')
    eje.grid(alpha=.2)
ax[1].legend(fontsize=7)
fig.tight_layout(); guardar_fig(fig,'comparacion_origen_comun'); plt.show(); plt.close(fig)

