"""Figuras trazables de perfiles, K y tau con theta fijo en 0.40."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
from datetime import datetime,timezone
from html import escape
import json,hashlib,sys,shutil,subprocess,os
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.windows import from_bounds
from rasterio.enums import Resampling
from scipy.integrate import cumulative_trapezoid
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
from scipy.stats import linregress
from pyproj import Transformer
from shapely.geometry import LineString,Point
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource,is_color_like,Normalize
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

CATEGORIAS={
 '01_Perfiles_chi_longitudinal':'Perfiles longitudinales y chi con knickpoints litologia y fallas',
 '02_ksn':'Perfiles de ksn de los 18 cauces',
 '03_Metodo_pendiente':'Calibracion K por pendiente y U de terraza',
 '04_PT_residuos_y_B_planta':'Interceptos PT residuos y seleccion manual de B en planta',
 '05_Metodo_George':'Calibracion K por el metodo de propagacion',
 '06_Tau_knickpoints':'Perfiles altura tau y tiempos de los knickpoints',
 '07_Dispersion_George':'Dispersion de K por propagacion',
 '08_Dispersion_pendiente':'Dispersion de K por pendiente',
 '09_Comparacion_metodos':'Comparacion de K por propagacion y pendiente',
 '10_Norte_sur_K_U':'Variacion norte sur de K y U por metodo',
 '11_Reloca_cuatro_cortes':'Reloca con cuatro cortes alternativos',
 '12_Mapa_regional_red':'Red regional con area aportante minima de 0.5 km2',
 '13_Tablas_K_tau':'Tablas de K promedio y tau de los knickpoints'}

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
 return h.hexdigest()

def guardar_json(p,data):Path(p).write_text(json.dumps(data,ensure_ascii=False,indent=2,default=str),encoding='utf-8')

def cargar(P):
 assert float(P['THETA'])==.40 and float(P['A0_M2'])==1e6 and int(P['N'])==1,'Esta entrega requiere theta=0.40, A0=1e6, n=1.'
 base=Path(P['BASE']);apoyo=base/'Codigo_consolidado_20261019';sys.path.insert(0,str(apoyo))
 import red_fallas as rf
 import revision_b_base as rb
 import calculo_unificado as cu
 from base_dos_K import project
 cfg=json.loads((base/'Configuracion.json').read_text(encoding='utf-8-sig'))
 out=Path(P['SALIDA']);datos=Path(P['DATOS']);out.mkdir(exist_ok=True);datos.mkdir(exist_ok=True)
 for name in CATEGORIAS:(out/name).mkdir(exist_ok=True)
 run=Path(P['CORRIDA_K']);k=pd.read_csv(run/'K_tau_U_todos_los_escenarios.csv',dtype={'mis':str},float_precision='round_trip')
 assert len(k)==24 and k.rol_resultado.eq('oficial').all()
 assert k[k.mis.eq('5c')].T_ka.eq(106).all() and k[k.mis.eq('5e')].T_ka.eq(125).all()
 names=pd.read_csv(base/'Productos/03_K_U/v2/resultados_K_U_v1.csv').set_index('rid').rio.to_dict()
 profiles={rid:pd.read_csv(run/f'perfil_R{rid:02d}_utilizado.csv',float_precision='round_trip') for rid in range(1,19)}
 chierr={}
 for rid,p in profiles.items():
  calc=1e6**.4*cumulative_trapezoid(p.area_m2.to_numpy()**-.4,p.distancia_m,initial=0)
  chierr[rid]=float(np.max(np.abs(calc-(p.chi_m-p.chi_m.iloc[0]))))
  if chierr[rid]>1e-6:raise ValueError(f'Chi del rio {rid} no corresponde a theta 0.40: {chierr[rid]}')
 geolpath=base/cfg['geologia_integrada'];geohash=sha(geolpath)
 cached=datos/'Perfiles_litologia';cached.mkdir(exist_ok=True)
 cachekey=datos/'Cache_litologia.json';expected={'geologia_sha':geohash,'perfiles':{str(r):sha(run/f'perfil_R{r:02d}_utilizado.csv') for r in profiles}}
 valid=cachekey.exists() and json.loads(cachekey.read_text(encoding='utf-8'))==expected and all((cached/f'R{r:02d}.csv').exists() for r in profiles)
 print('Cargando litología de 18 perfiles'+(' desde caché.' if valid else ' por intersección cartográfica.'),flush=True)
 if valid:
  profiles={r:pd.read_csv(cached/f'R{r:02d}.csv',float_precision='round_trip').fillna({'codigo_lito':'Sin cobertura','color_lito':'#8b949e','descripcion_lito':'Sin cobertura cartográfica'}) for r in profiles}
  geology=gpd.read_file(geolpath).to_crs(32718)
 else:
  profiles,geology=rf.asignar_litologia(profiles,geolpath)
  for rid,p in profiles.items():p.to_csv(cached/f'R{int(rid):02d}.csv',index=False,encoding='utf-8-sig')
  guardar_json(cachekey,expected)
 profiles={int(r):p for r,p in profiles.items()}
 for p in profiles.values():
  p['color_lito']=[P['COLORES_LITO'].get(str(c),v if is_color_like(v) else '#8b949e') for c,v in zip(p.codigo_lito,p.color_lito)]
  stats=rf.ventana_ols(p.distancia_m,p.chi_m,p.elevacion_m,W=P['VENTANA_KSN_M'],min_n=P['MIN_NODOS_KSN'],min_fraction=P['FRACCION_VENTANA_MIN'],theta=.4,A0=1e6)
  p['ksn_040_figuras']=stats['ksn'];p['R2_ksn_figuras']=stats['R2'];p['ventana_efectiva_m']=stats['span']
 faults=gpd.read_file(base/cfg['chaf']).to_crs(32718);faults['f_id']=pd.to_numeric(faults.f_id).astype(int)
 faults['nombre_base']=faults.falla.astype(str).str.replace('*','',regex=False).str.strip()
 crossings=rf.cruces_en_perfiles(profiles,faults)
 if crossings.empty:crossings=pd.DataFrame(columns=['perfil_id','f_id','falla','nombre_base','chi_cruce_m','dist_cruce_m'])
 else:
  crossings['qx']=(crossings.x_cruce*100).round();crossings['qy']=(crossings.y_cruce*100).round()
  crossings=crossings.drop_duplicates(['perfil_id','nombre_base','qx','qy']).drop(columns=['qx','qy'])
 manual=pd.read_csv(P.get('CSV_KP_ORIGINAL',base/cfg['kp_manual_csv']),sep=cfg['kp_separador'],decimal=cfg['kp_decimal']);kps=[]
 for row in manual.itertuples():
  if not np.isfinite([row.id,row.x_utm,row.y_utm]).all():continue
  rid=int(row.id)
  if rid not in profiles:continue
  at=project(profiles[rid],[row.x_utm,row.y_utm])
  kps.append({'rid':rid,'id_kp':row.id_kp,'theta_seleccion_original':row.theta_seleccion,**at})
 kp=pd.DataFrame(kps)
 # Recalcular ambos métodos desde los mismos datos para verificar las tablas vigentes.
 fitrows=[]
 for v in k.itertuples():
  p=profiles[int(v.rid)];ss,xx,zz,n=cu.tramo_ajuste(p,v.s_B_m,v.s_PT_m);fit=linregress(xx,zz)
  kg=(v.chi_PT_m-v.chi_B_m)/(1e6**.4*v.T_ka*1000)
  ksn=fit.slope*1e6**.4;up=(v.H_publicada_para_U_m-v.nivel_mar_para_U_m)/v.T_ka
  ks=(up/1000)/ksn
  assert np.isclose(kg,v.K_m02_anio,rtol=1e-12) and np.isclose(ks,v.K_pendiente_terraza_m02_anio,rtol=1e-12)
  fitrows.append({'id_medicion':v.id_medicion,'M_chi':fit.slope,'intercepto_z_m':fit.intercept,'R2':fit.rvalue**2,
     'RMSE_m':float(np.sqrt(np.mean((zz-fit.slope*xx-fit.intercept)**2))),
     'K_George':kg,'K_pendiente':ks,'ksn_ajuste':ksn,'U_terraza':up,'n_nodos_reales':n})
 fits=pd.DataFrame(fitrows);k=k.merge(fits,on='id_medicion',validate='one_to_one')
 k['U_George']=k.K_George*k.ksn_ajuste*1000
 k['U_pendiente']=k.K_pendiente*k.ksn_ajuste*1000
 k['tau_KP_pendiente_ka']=np.nan
 for idx,v in k.iterrows():
  rr=kp[kp.rid.eq(int(v.rid))]
  if len(rr):
   last=rr.loc[rr.s_m.idxmax()]
   k.loc[idx,'tau_KP_pendiente_ka']=(last.chi_m-v.chi_B_m)/(v.K_pendiente*1e6**.4*1000)
 # Para promedios descriptivos no se pondera tres veces un mismo P_T físico.
 unique=k.drop_duplicates(['rid','mis','s_PT_m','chi_PT_m'])
 ns=unique.groupby(['rid','rio','mis']).agg(latitud_sur=('latitud_sur','first'),
    K_George=('K_George','mean'),K_pendiente=('K_pendiente','mean'),U_George=('U_George','mean'),U_pendiente=('U_pendiente','mean'),
    n_PT_distintos=('id_medicion','count')).reset_index()
 # Exportar tablas científicas que sostienen las figuras, separadas de la apariencia.
 k.to_csv(datos/'Calibraciones_figuras.csv',index=False,encoding='utf-8-sig');kp.to_csv(datos/'KP_proyectados_theta040.csv',index=False,encoding='utf-8-sig')
 crossings.to_csv(datos/'Cruces_fallas_18_perfiles.csv',index=False,encoding='utf-8-sig');ns.to_csv(datos/'Promedios_N_S_por_rio_MIS.csv',index=False,encoding='utf-8-sig')
 for rid,p in profiles.items():p.to_csv(datos/f'R{rid:02d}_perfil_figuras.csv',index=False,encoding='utf-8-sig')
 C={'P':P,'base':base,'out':out,'datos':datos,'run':run,'profiles':profiles,'names':names,'k':k,'kp':kp,
    'geology':geology,'faults':faults,'crossings':crossings,'ns':ns,'dem':base/cfg['dem'],'cfg':cfg}
 trace={'fecha_utc':datetime.now(timezone.utc).isoformat(),'theta':.4,'A0':1e6,'n':1,'n_rios_perfiles':18,'n_rios_calibrados':6,
    'n_casos_B_MIS':8,'n_escenarios':24,'n_KP':len(kp),'max_error_chi_theta040_m':max(chierr.values()),
    'corrida_K':str(run),'fuentes_sha256':{str(f):sha(f) for f in [run/'K_tau_U_todos_los_escenarios.csv',geolpath,base/cfg['chaf'],Path(P.get('CSV_KP_ORIGINAL',base/cfg['kp_manual_csv']))]},
    'metodo_pendiente':'K=(U_terraza/1000)/(Mchi*A0**theta); U_terraza=(H_catalogo-nivel_mar)/T_publicado',
    'metodo_George':'K=delta_chi/(A0**theta*T_publicado*1000)',
    'promedios_NS':'Media de PT físicos distintos por río y MIS, sin mezclar 5c y 5e',
    'U':'George: U aparente=1000*K*ksn; pendiente: U de terraza impuesto, no estimación independiente',
    'parametros':P}
 guardar_json(datos/'Trazabilidad_figuras.json',trace);guardar_json(out/'Parametros_figuras.json',P)
 plt.rcParams.update({'font.family':P['FUENTE'],'font.size':P['TAMANO_FUENTE'],'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 print(f'Listos: 18 perfiles; {len(kp)} KP manuales; 8 casos oficiales; theta=0.40.',flush=True)
 return C

def iniciar(C,cat):
 f=C['out']/'Catalogo_figuras.json'
 entries=json.loads(f.read_text(encoding='utf-8')) if f.exists() else []
 entries=[e for e in entries if e['categoria']!=cat];guardar_json(f,entries)

def guardar(C,fig,cat,name,title,caption,rid=None,mis=None):
 if cat in ['03_Metodo_pendiente','04_PT_residuos_y_B_planta','05_Metodo_George'] and rid==9:caption+=' R09: se conserva el perfil y K de la calibración original; la cuenca recortada P5 se aplica a los perfiles vigentes y tiempos, sin recalibrar K.'
 stem=C['out']/cat/name
 if C['P']['ESCALA_FIGURAS']!=1:fig.set_size_inches(np.asarray(fig.get_size_inches())*C['P']['ESCALA_FIGURAS'])
 for ext in C['P']['FORMATOS']:fig.savefig(stem.with_suffix('.'+ext),dpi=C['P']['DPI'],facecolor='white',bbox_inches='tight')
 plt.close(fig)
 f=C['out']/'Catalogo_figuras.json';entries=json.loads(f.read_text(encoding='utf-8')) if f.exists() else []
 entries.append({'categoria':cat,'nombre':name,'titulo':title,'leyenda':caption,'rid':rid,'mis':mis,
    'png':str(stem.with_suffix('.png').relative_to(C['out']).as_posix()),'formatos':C['P']['FORMATOS'],
    'parametros_sha256':hashlib.sha256(json.dumps(C['P'],sort_keys=True).encode()).hexdigest()})
 guardar_json(f,entries)

def color(C,mis):return C['P']['COLORES_MIS'][str(mis)]

def suavizar(C,y,s):
 a=np.asarray(y,float).copy()
 if not C['P']['SUAVIZAR_PERFILES']:return a
 good=np.isfinite(a);ix=np.flatnonzero(good)
 for chunk in np.split(ix,np.flatnonzero(np.diff(ix)>1)+1):
  if len(chunk)<7:continue
  step=np.median(np.diff(np.asarray(s)[chunk]));w=int(C['P']['VENTANA_SUAVIZADO_M']/max(step,1e-6));w=max(5,w//2*2+1);w=min(w,len(chunk)//2*2-1)
  if w>=5:a[chunk]=savgol_filter(a[chunk],w,2)
 return a

def titulo(C,rid):return f'R{int(rid):02d} · {C["names"][int(rid)]}'

def cruces(C,rid,ax,xmode='chi'):
 if not C['P']['MOSTRAR_FALLAS']:return
 q=C['crossings'][C['crossings'].perfil_id.eq(rid)]
 for j,v in enumerate(q.itertuples()):
  x=v.chi_cruce_m if xmode=='chi' else v.dist_cruce_m/1000
  co=C['P']['COLORES_FALLAS'].get(v.nombre_base,'#444444');ax.axvline(x,ls=':',c=co,lw=C['P']['ANCHO_FALLA'],alpha=.9)
  ax.text(x,.98-(j%3)*.15,v.nombre_base,transform=ax.get_xaxis_transform(),rotation=90,va='top',ha='right',fontsize=7,color=co)

def marcar_KP(C,rid,ax,xmode='chi'):
 q=C['kp'][C['kp'].rid.eq(rid)]
 for v in q.itertuples():
  x=v.chi_m if xmode=='chi' else v.s_m/1000
  ax.scatter(x,v.z_m,c=C['P']['COLOR_KP'],s=C['P']['TAMANO_KP'],marker=C['P']['SIMBOLOS']['KP'],zorder=6,edgecolors='white',linewidths=.5)
  ax.annotate('KP',(x,v.z_m),xytext=(5,8),textcoords='offset points',fontsize=8,color=C['P']['COLOR_KP'])

def barra_lito(C,p,ax,x):
 if not C['P']['MOSTRAR_LITOLOGIA']:ax.set_axis_off();return []
 codes=p.codigo_lito.to_numpy();bound=np.r_[0,np.flatnonzero(codes[1:]!=codes[:-1])+1,len(p)-1]
 for start,end in zip(bound[:-1],bound[1:]):
  right=float(x[end]);left=float(x[start]);ax.axvspan(left,right,color=p.color_lito.iloc[start],lw=0)
 ax.set_ylim(0,1);ax.set_xlim(float(x[0]),float(x[-1]));ax.set_yticks([]);ax.tick_params(axis='x',bottom=False,labelbottom=False)
 ax.set_ylabel('Lito.',rotation=0,ha='right',va='center',fontsize=8)
 return [Patch(fc=g.color_lito.iloc[0],label=str(code)) for code,g in p.groupby('codigo_lito',sort=False)]

def perfiles(C):
 cat='01_Perfiles_chi_longitudinal';iniciar(C,cat)
 for rid,p in C['profiles'].items():
  fig=plt.figure(figsize=(13,6.5));gs=fig.add_gridspec(2,2,height_ratios=[10,.42],hspace=.35,wspace=.23)
  ax1=fig.add_subplot(gs[0,0]);ax2=fig.add_subplot(gs[0,1]);b1=fig.add_subplot(gs[1,0]);b2=fig.add_subplot(gs[1,1])
  z=suavizar(C,p.elevacion_m,p.distancia_m)
  for ax,x,mode,label in [(ax1,p.distancia_m.to_numpy()/1000,'dist','Distancia aguas arriba (km)'),(ax2,p.chi_m.to_numpy(),'chi','χ (m)')]:
   ax.plot(x,z,c=C['P']['COLOR_PERFIL'],lw=1.25);marcar_KP(C,rid,ax,mode);cruces(C,rid,ax,mode)
   ax.set_xlabel(label);ax.set_ylabel('Altura DEM (m)');ax.grid(alpha=.18);ax.set_xlim(float(x[0]),float(x[-1]))
  handles=barra_lito(C,p,b1,p.distancia_m.to_numpy()/1000);barra_lito(C,p,b2,p.chi_m.to_numpy())
  if handles:fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.005),ncol=min(8,len(handles)),fontsize=8)
  fig.subplots_adjust(bottom=.19,top=.88);fig.suptitle(titulo(C,rid)+' · perfil longitudinal y χ · θ=0,40',y=.98)
  missing=float(p.codigo_lito.eq('Sin cobertura').mean());nkp=int(C['kp'].rid.eq(rid).sum())
  sub=f'{nkp} KP manual(es) · sin cobertura litológica: {missing:.0%}'
  if not nkp:sub+=' · no se asigna un KP automático'
  fig.text(.5,.91,sub,ha='center',fontsize=9)
  guardar(C,fig,cat,f'R{rid:02d}_longitudinal_chi',titulo(C,rid)+' Perfiles longitudinal y chi',
    'Perfil del DEM original y χ con θ=0,40. Estrellas: KP manuales. Franjas: unidades cartográficas; la etiqueta Sin cobertura identifica falta de cartografía. Los colores de las otras unidades se reconocen por su código en la leyenda. Líneas punteadas: intersecciones reales con trazas CHAF. Una coincidencia espacial no demuestra control litológico o actividad tectónica.',rid)

def ksn(C):
 cat='02_ksn';iniciar(C,cat)
 for rid,p in C['profiles'].items():
  fig=plt.figure(figsize=(11,5.8));gs=fig.add_gridspec(2,1,height_ratios=[10,.5],hspace=.35)
  ax=fig.add_subplot(gs[0]);band=fig.add_subplot(gs[1]);x=p.distancia_m.to_numpy()/1000
  ax.plot(x,suavizar(C,p.ksn_040_figuras,p.distancia_m),c=C['P']['COLOR_KSN'],lw=1.1)
  for v in C['kp'][C['kp'].rid.eq(rid)].itertuples():ax.axvline(v.s_m/1000,c=C['P']['COLOR_KP'],ls='--',lw=1,label='KP manual')
  cruces(C,rid,ax,'dist');handles=barra_lito(C,p,band,x)
  ax.set_xlabel('Distancia aguas arriba (km)');ax.set_ylabel(r'$k_{sn}$ (m$^{0,8}$)');ax.grid(alpha=.18);ax.set_xlim(x[0],x[-1])
  if C['P']['KSN_LOG']:ax.set_yscale('log')
  if handles:fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.005),ncol=min(8,len(handles)),fontsize=8)
  fig.subplots_adjust(bottom=.20,top=.85);fig.suptitle(titulo(C,rid)+f' · ksn local · ventana {C["P"]["VENTANA_KSN_M"]:g} m · θ=0,40')
  guardar(C,fig,cat,f'R{rid:02d}_ksn',titulo(C,rid)+' ksn local',
    f'OLS altura–χ con intercepto libre en ventana móvil de {C["P"]["VENTANA_KSN_M"]:g} m. ksn=Mχ A₀^θ; θ=0,40. Pendientes no positivas, ventanas insuficientes y huecos permanecen inválidos. El suavizado opcional sólo afecta la apariencia, no los ajustes de K.',rid)

def tramo(C,v):
 import calculo_unificado as cu
 p=C['profiles'][int(v.rid)];return cu.tramo_ajuste(p,v.s_B_m,v.s_PT_m)[:3]

def metodo_pendiente(C):
 cat='03_Metodo_pendiente';iniciar(C,cat)
 for case,q in C['k'].groupby('caso_B',sort=False):
  q=q.sort_values('orden_PT');v=q[q.escenario_PT.eq('intermedio_real')].iloc[0];p=C['profiles'][int(v.rid)]
  fig,(ax,ar)=plt.subplots(2,1,figsize=(10.5,7),gridspec_kw={'height_ratios':[3,1]},layout='constrained')
  mask=p.distancia_m.between(v.s_B_m,float(q.s_PT_m.max())+100)
  ax.plot(p.loc[mask,'chi_m'],p.loc[mask,'elevacion_m'],c='#777',lw=1,label='Perfil DEM')
  for j,row in enumerate(q.itertuples()):
   ss,xx,zz=tramo(C,row);co=color(C,row.mis);style=[':', '-', '--'][j]
   fit=row.M_chi*xx+row.intercepto_z_m
   ax.plot(xx,fit,c=co,ls=style,lw=1.35,label=f'{row.escenario_PT}: K={row.K_pendiente:.3e}')
   ax.scatter(row.chi_PT_m,row.H_m,c=co,marker=C['P']['SIMBOLOS']['PT'],s=40)
   ar.plot(xx,zz-fit,c=co,ls=style,lw=.9)
  ax.set_ylabel('Altura DEM (m)');ar.set_ylabel('Residuo OLS (m)');ar.set_xlabel('χ (m)');ar.axhline(0,c='#777',lw=.8)
  ax.legend(loc='lower right',fontsize=8);ax.grid(alpha=.2);ar.grid(alpha=.2)
  note=(f'Mχ intermedio = {v.M_chi:.5f} · R² = {v.R2:.3f}\n'
        f'ksn = Mχ A₀^θ = {v.ksn_ajuste:.2f} m⁰·⁸\n'
        f'U terraza = (H catálogo − nivel mar)/T = {v.U_terraza:.3f} mm/año\n'
        f'K pendiente = (U terraza/1000)/ksn\nT = {v.T_ka:g} ± {v.error_T_ka:g} ka · θ = 0,40')
  ax.text(.025,.98,note,transform=ax.transAxes,va='top',fontsize=9,bbox={'fc':'white','ec':'none','alpha':.90})
  fig.suptitle(titulo(C,v.rid)+f' · {v.id_B_corto} · MIS {v.mis} · método de pendientes')
  guardar(C,fig,cat,f'R{int(v.rid):02d}_{v.id_B_corto}_{v.mis}_K_pendiente',titulo(C,v.rid)+f' Pendiente MIS {v.mis}',
    'Ajuste OLS altura–χ con intercepto libre para cada tramo B–Pₜ. K=(U_terraza/1000)/(Mχ A₀^θ), con U de terraza calculado usando altura del catálogo, nivel marino de referencia y edad publicada. El panel inferior muestra residuos del ajuste. Los tres interceptos se conservan por separado.',int(v.rid),v.mis)

def metodo_george(C):
 cat='05_Metodo_George';iniciar(C,cat)
 for case,q in C['k'].groupby('caso_B',sort=False):
  q=q.sort_values('orden_PT');v=q.iloc[0];p=C['profiles'][int(v.rid)]
  mask=p.distancia_m.between(v.s_B_m,float(q.s_PT_m.max())+150)
  fig,ax=plt.subplots(figsize=(10.5,6),layout='constrained')
  ax.plot(p.loc[mask,'chi_m']-v.chi_B_m,p.loc[mask,'elevacion_m'],c=C['P']['COLOR_PERFIL'],lw=1.4)
  ax.axhline(v.H_m,c=color(C,v.mis),ls=':',lw=1,label=f'H DEM terraza = {v.H_m:.2f} m')
  ax.scatter(0,v.z_B_m,c='#222',marker=C['P']['SIMBOLOS']['B'],s=60,zorder=5,label=f'B manual {v.id_B_corto}')
  for j,row in enumerate(q.itertuples()):
   dx=row.chi_PT_m-row.chi_B_m;ax.axvline(dx,color=color(C,row.mis),ls=[':', '-', '--'][j],alpha=.55,lw=.8)
   ax.scatter(dx,row.H_m,c=color(C,row.mis),marker=C['P']['SIMBOLOS']['PT'],s=45,zorder=5)
  labels=['K George = [χ(Pₜ) − χ(B)] / [A₀^θ T]','',f'T = {v.T_ka:g} ± {v.error_T_ka:g} ka · θ = 0,40']
  for row in q.itertuples():labels.append(f'{row.escenario_PT}: Δχ={row.delta_chi_m:.2f} m; K={row.K_George:.4e}')
  ax.text(.025,.98,'\n'.join(labels),transform=ax.transAxes,va='top',fontsize=9,bbox={'fc':'white','ec':'none','alpha':.9})
  ax.set_xlabel('χ − χ(B) (m)');ax.set_ylabel('Altura DEM (m)');ax.grid(alpha=.2);ax.legend(loc='lower right',fontsize=8)
  fig.suptitle(titulo(C,v.rid)+f' · {v.id_B_corto} · MIS {v.mis} · método de propagación')
  guardar(C,fig,cat,f'R{int(v.rid):02d}_{v.id_B_corto}_{v.mis}_K_George',titulo(C,v.rid)+f' George MIS {v.mis}',
    'B se fija manualmente y H corresponde al DEM en la referencia de terraza. Los cruces reales H–perfil definen Pₜ primero, intermedio real y último. Cada diferencia χ(Pₜ)−χ(B) y edad T produce un K; las líneas verticales ubican los interceptos. No se ajusta K promediando MIS.',int(v.rid),v.mis)

def relieve(C,bounds,maxside=700):
 with rasterio.open(C['dem']) as src:
  b=(max(bounds[0],src.bounds.left),max(bounds[1],src.bounds.bottom),min(bounds[2],src.bounds.right),min(bounds[3],src.bounds.top))
  win=from_bounds(*b,src.transform);w=max(2,int(win.width));h=max(2,int(win.height));scale=min(1,maxside/max(w,h));shape=(max(2,int(h*scale)),max(2,int(w*scale)))
  z=src.read(1,window=win,out_shape=shape,masked=True,resampling=Resampling.bilinear).astype(float).filled(np.nan)
 good=np.isfinite(z);fill=np.where(good,z,np.nanmedian(z));dx=(b[2]-b[0])/shape[1];dy=(b[3]-b[1])/shape[0]
 hs=LightSource(azdeg=315,altdeg=45).hillshade(fill,vert_exag=1,dx=dx,dy=dy);hs[~good]=np.nan
 return z,hs,b

def residuos_y_planta(C):
 cat='04_PT_residuos_y_B_planta';iniciar(C,cat)
 catalog=pd.read_csv(C['base']/'Codigo_consolidado_20261019/Fuentes/terrazas_170_UTM18S_con_edad_referencia_Jara_Munoz_2015.csv')
 for case,q in C['k'].groupby('caso_B',sort=False):
  q=q.sort_values('orden_PT');v=q.iloc[0];p=C['profiles'][int(v.rid)];pad=C['P']['RADIO_PLANTA_B_M']
  bounds=(v.x_B-pad,v.y_B-pad,v.x_B+pad,v.y_B+pad);z,hs,b=relieve(C,bounds)
  fig,axes=plt.subplots(1,3,figsize=(16,5.5),layout='constrained');am,ap,ar=axes
  am.imshow(hs,extent=np.array(b)[[0,2,1,3]]/1000,origin='upper',cmap='gray',vmin=0,vmax=1,alpha=.8)
  local=p[p.x_utm.between(b[0],b[2])&p.y_utm.between(b[1],b[3])]
  am.plot(local.x_utm/1000,local.y_utm/1000,c='#2081b3',lw=1.4,label='Cauce actual')
  for row in C['k'][C['k'].rid.eq(v.rid)].drop_duplicates('id_B_corto').itertuples():
   am.scatter(row.x_B/1000,row.y_B/1000,c=color(C,row.mis),s=70,marker=C['P']['SIMBOLOS']['B'],edgecolors='white',zorder=6,label=f'B {row.id_B_corto} · {row.mis}')
  nearby=catalog[catalog.este_utm18s_m.between(b[0],b[2])&catalog.norte_utm18s_m.between(b[1],b[3])]
  for mis,symbol in [('5e','^'),('5c','v')]:
   a=nearby[nearby.mis_asignado.astype(str).str.lower().str.replace(' ','').str.contains(mis)]
   if len(a):am.scatter(a.este_utm18s_m/1000,a.norte_utm18s_m/1000,c=color(C,mis),s=35,marker=symbol,label=f'Terrazas {mis}')
  am.scatter(v.x_terraza/1000,v.y_terraza/1000,facecolors='none',edgecolors='#222',s=100,zorder=7,label='Referencia usada')
  if C['P']['MOSTRAR_FALLAS']:C['faults'].plot(ax=am,color='#555',linestyle=':',linewidth=.8)
  am.set_xlim(b[0]/1000,b[2]/1000);am.set_ylim(b[1]/1000,b[3]/1000);am.set_aspect('equal');am.set_xlabel('Este UTM 18S (km)');am.set_ylabel('Norte UTM 18S (km)');am.set_title('Ubicación del B seleccionado manualmente');am.legend(fontsize=6.8,loc='best')
  am.annotate('N',xy=(.94,.94),xytext=(.94,.80),xycoords='axes fraction',ha='center',arrowprops={'arrowstyle':'->','color':'black'})
  # Los huecos de selección de puntos no se unen artificialmente en los datos.
  show=p[p.distancia_m.between(v.s_B_m,float(q.s_PT_m.max())+250)]
  ap.plot(show.chi_m-v.chi_B_m,show.elevacion_m,c='#536f80',lw=1.1);ap.axhline(v.H_m,c=color(C,v.mis),ls='--',lw=.9)
  ap.scatter(0,v.z_B_m,c='#222',marker=C['P']['SIMBOLOS']['B'],s=45)
  ptlo,ptmax=float(q.s_PT_m.min()),float(q.s_PT_m.max());span=max(350,ptmax-ptlo)
  sub=p[p.distancia_m.between(ptlo-.25*span,ptmax+.25*span)]
  ar.plot((sub.distancia_m-ptlo)/1000,sub.elevacion_m-v.H_m,c='#536f80',lw=1.1);ar.axhline(0,c='#222',lw=.8)
  for j,row in enumerate(q.itertuples()):
   ap.scatter(row.chi_PT_m-row.chi_B_m,row.H_m,c=color(C,v.mis),s=45,marker=C['P']['SIMBOLOS']['PT'])
   ar.scatter((row.s_PT_m-ptlo)/1000,0,c=color(C,v.mis),s=40,marker=['o','D','s'][j],label=row.escenario_PT)
  ap.set_xlabel('χ − χ(B) (m)');ap.set_ylabel('Altura DEM (m)');ap.set_title(f'H = {v.H_m:.2f} m e interceptos Pₜ')
  ar.set_xlabel('Distancia desde el primer Pₜ (km)');ar.set_ylabel('Z perfil − H (m)');ar.set_title('Residuo de intersección con H');ar.legend(fontsize=8)
  ap.grid(alpha=.2);ar.grid(alpha=.2);fig.suptitle(titulo(C,v.rid)+f' · B={v.id_B_corto} · MIS {v.mis}')
  guardar(C,fig,cat,f'R{int(v.rid):02d}_{v.id_B_corto}_{v.mis}_B_PT_residuo',titulo(C,v.rid)+f' B y PT MIS {v.mis}',
    'Izquierda: posición manual de B en el relieve y referencias de terraza por MIS; B no se determina automáticamente mediante una isohipsa. Centro: H del DEM y cruces Pₜ. Derecha: residuo Z−H, cuya raíz define cada cruce. Los cruces exactos tienen residuo cero por interpolación, lo que no valida su equivalencia cronológica.',int(v.rid),v.mis)

def tau(C):
 cat='06_Tau_knickpoints';iniciar(C,cat);rows=[]
 for rid,g in C['k'].groupby('rid'):
  p=C['profiles'][int(rid)];kp=C['kp'][C['kp'].rid.eq(rid)]
  maxs=float(kp.s_m.max()) if len(kp) else float(p.distancia_m.iloc[-1])
  if not C['P']['TAU_HASTA_KP']:maxs=float(p.distancia_m.iloc[-1])
  fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True,layout='constrained')
  for ax,method in zip(axes,['George','pendiente']):
   for case,q in g.groupby('caso_B',sort=False):
    v=q.iloc[0];s=np.unique(np.r_[v.s_B_m,p.loc[p.distancia_m.between(v.s_B_m,maxs),'distancia_m'],maxs]);xx=np.interp(s,p.distancia_m,p.chi_m);zz=np.interp(s,p.distancia_m,p.elevacion_m)
    kval=q['K_'+method].to_numpy();matrix=np.array([(xx-v.chi_B_m)/(kv*1e6**.4*1000) for kv in kval])
    for j,row in enumerate(q.sort_values('orden_PT').itertuples()):
     kval=getattr(row,'K_'+method);tt=(xx-row.chi_B_m)/(kval*1e6**.4*1000)
     ax.plot(tt,zz,c=color(C,row.mis),ls=[':', '-', '--'][j],lw=1,label=f'{row.mis} · {row.id_B_corto} · {row.escenario_PT}')
     for point in kp.itertuples():
      if point.chi_m<row.chi_B_m:continue
      tk=(point.chi_m-row.chi_B_m)/(kval*1e6**.4*1000)
      ax.scatter(tk,point.z_m,c=color(C,row.mis),marker=C['P']['SIMBOLOS']['KP'],s=55,zorder=5)
      rows.append({'rid':int(rid),'rio':row.rio,'mis':row.mis,'id_medicion':row.id_medicion,'metodo':method,'id_kp':point.id_kp,'tau_KP_ka':tk})
   ax.set_xlabel('τ desde B de cada escenario (ka)');ax.grid(alpha=.2);ax.set_title('K por '+method);ax.set_xlim(left=0);ax.legend(fontsize=6.8,loc='lower right')
  axes[0].set_ylabel('Altura DEM (m)');fig.suptitle(titulo(C,rid)+' · altura–τ · K uniforme por escenario')
  if not len(kp):fig.suptitle(titulo(C,rid)+' · altura–τ · K uniforme por escenario\nSin KP manual registrado: se muestra el perfil completo sin asignar una edad a un KP.')
  guardar(C,fig,cat,f'R{int(rid):02d}_altura_tau',titulo(C,rid)+' Tau de knickpoints',
    'τ=(χ−χB)/(K A₀^θ). Cada línea usa un K calibrado y su B correspondiente; las estrellas ubican los KP manuales. Los paneles separan los dos métodos de K. Las líneas distinguen las elecciones de Pₜ a edad nominal, sin áreas rellenas. Son tiempos de propagación condicionales; no dataciones tectónicas directas. En Chanco no se inventa un KP ausente del registro.',int(rid))
 pd.DataFrame(rows).to_csv(C['datos']/'Tau_KP_dos_metodos.csv',index=False,encoding='utf-8-sig')

def dispersion(C,metodos):
 cat={'George':'07_Dispersion_George','pendiente':'08_Dispersion_pendiente','ambos':'09_Comparacion_metodos'}[metodos];iniciar(C,cat)
 groups=[(int(r),g) for r,g in C['k'].groupby('rid')]+[(None,C['k'])]
 for rid,q in groups:
  q=q.sort_values(['rid','mis','orden_PT']);fig,ax=plt.subplots(figsize=(max(10,.4*len(q)),5.8),layout='constrained')
  xpos={v:i for i,v in enumerate(q.id_medicion)};methods=['George','pendiente'] if metodos=='ambos' else [metodos]
  for m,method in enumerate(methods):
   for mis,g in q.groupby('mis'):
    offset=(m-.5)*.22 if len(methods)==2 else 0
    ax.scatter([xpos[x]+offset for x in g.id_medicion],g['K_'+method]*1e5,c=color(C,mis),s=48,
      marker=C['P']['SIMBOLOS'][method],label=f'{method} · MIS {mis}',alpha=.9,zorder=4)
    distinct=g.drop_duplicates(['rid','mis','s_PT_m','chi_PT_m'])
    mean=float(distinct['K_'+method].mean())*1e5
    ax.axhline(mean,c=color(C,mis),ls='--' if method=='George' else '-.',lw=1.2,
      label=f'Media {method} · {mis}: {mean:.3f} ×10⁻⁵')
  ax.set_xticks(range(len(q)),q.id_medicion,rotation=65,ha='right',fontsize=8 if len(q)<10 else 7)
  ax.set_ylabel('K (×10⁻⁵ m⁰·²/año)');ax.set_xlabel('ID de escenario · orden sin significado espacial');ax.grid(axis='y',alpha=.2);ax.legend(fontsize=9)
  title=(titulo(C,rid) if rid is not None else 'Seis cauces oficiales')+' · K por '+('ambos métodos' if metodos=='ambos' else metodos)
  fig.suptitle(title);name=f'R{rid:02d}' if rid is not None else 'Todos_los_rios'
  guardar(C,fig,cat,name+'_dispersion_K',title,
    'Cada ID identifica río, B, MIS y etiqueta Pₜ. Los colores separan MIS 5c y 5e; las formas separan métodos cuando se comparan. Los tres Pₜ coincidentes se conservan como etiquetas y no como observaciones independientes. El eje X sólo sirve para inspeccionar la dispersión. Las líneas horizontales muestran la media de Pₜ físicos distintos dentro de cada MIS y método; no son ajustes ni mezclan MIS.',rid)

def norte_sur(C):
 cat='10_Norte_sur_K_U';iniciar(C,cat)
 trendrows=[]
 for method in ['George','pendiente']:
  fig,axes=plt.subplots(2,2,figsize=(12.5,9),sharex='col',layout='constrained')
  for col,mis in enumerate(['5e','5c']):
   q=C['ns'][C['ns'].mis.eq(mis)].sort_values('latitud_sur')
   x=q.latitud_sur.to_numpy()
   for row,var in enumerate(['K','U']):
    ax=axes[row,col];y=q[var+'_'+method].to_numpy();scale=1e5 if var=='K' else 1
    ax.scatter(x,y*scale,s=50,c=color(C,mis),marker=C['P']['SIMBOLOS'][method],zorder=4)
    if len(x)>1:ax.plot(x,y*scale,c=color(C,mis),lw=1.1)
    trendrows.extend({'metodo':method,'mis':mis,'variable':var,'latitud_sur':float(a),'valor_punto':float(b)} for a,b in zip(x,y))
    for j,(a,b,label) in enumerate(zip(x,y*scale,q.rio)):ax.annotate(label.replace('Estero de ','').replace('Estero ','').replace('Rio ',''),(a,b),xytext=(4,7 if j%2==0 else -13),textcoords='offset points',fontsize=7)
    ax.set_ylabel('K (×10⁻⁵ m⁰·²/año)' if var=='K' else ('U aparente (mm/año)' if method=='George' else 'U terraza impuesto (mm/año)'))
    ax.grid(alpha=.2);ax.margins(x=.13,y=.22)
    if row==0:ax.set_title(f'MIS {mis} · '+('seis ríos' if mis=='5e' else 'Chanco y Curanilahue')+' · media por río')
    else:ax.set_xlabel('Latitud (°S) · norte → sur')
  fig.suptitle('Variación norte–sur · método '+method+'\nPromedios descriptivos de Pₜ distintos, separados por MIS')
  guardar(C,fig,cat,'Norte_sur_K_U_'+method,'K y U de norte a sur por '+method,
    'Un promedio aritmético de los Pₜ físicos retenidos por río y MIS; no se mezclan 5c y 5e. Los segmentos unen directamente los promedios observados en orden de latitud, sin regresión ni suavizado. Para George, U=1000 K ksn es una tasa aparente; para pendiente, U se impone desde la terraza. Con dos ríos en 5c se dibuja una conexión lineal.')
 pd.DataFrame(trendrows).to_csv(C['datos']/'Guias_visuales_N_S.csv',index=False,encoding='utf-8-sig')

def regional(C):
 cat='12_Mapa_regional_red';iniciar(C,cat)
 src=Path(C['P']['RED_REGIONAL']);layer=C['P']['CAPA_RED'];poly=gpd.read_file(src,layer='poligono').to_crs(32718)
 amin=float(C['P']['AREA_MINIMA_RED_KM2']);assert amin>=.1,'La red fuente sólo contiene áreas desde 0.1 km².'
 edges=gpd.read_file(src,layer=layer,where=f'area_km2 >= {amin:.12g}').to_crs(32718)
 assert len(edges) and edges.area_km2.min()>=amin-1e-12
 # La geometría viene recortada al polígono; no se submuestrean las ramas.
 bounds=poly.total_bounds;pad=max(bounds[2]-bounds[0],bounds[3]-bounds[1])*.015
 b=(bounds[0]-pad,bounds[1]-pad,bounds[2]+pad,bounds[3]+pad)
 zz,hs,bb=relieve(C,b,maxside=1800)
 fig,ax=plt.subplots(figsize=(10.5,12),layout='constrained')
 if C['P']['MOSTRAR_RELIEVE_MAPA']:ax.imshow(hs,extent=np.array(bb)[[0,2,1,3]]/1000,origin='upper',cmap='gray',alpha=.6,vmin=0,vmax=1)
 from shapely.geometry import mapping
 segments=[];areas=[]
 for geom,a in zip(edges.geometry,edges.area_km2):
  for line in (geom.geoms if geom.geom_type=='MultiLineString' else [geom]):
   segments.append(np.array(line.coords)/1000);areas.append(a)
 widths=np.clip(.25+.28*np.log10(np.maximum(areas,amin)/amin+1),.25,1.5)
 collection=LineCollection(segments,colors=C['P']['COLOR_RED'],linewidths=widths,alpha=.8,rasterized=True,zorder=2)
 ax.add_collection(collection)
 polykm=poly.copy();polykm.geometry=polykm.geometry.scale(xfact=.001,yfact=.001,origin=(0,0));polykm.boundary.plot(ax=ax,color='#222',linewidth=1.2,zorder=4)
 if C['P']['MOSTRAR_FALLAS_MAPA']:
  fk=C['faults'].copy();fk.geometry=fk.geometry.scale(xfact=.001,yfact=.001,origin=(0,0))
  for name,g in fk.groupby('nombre_base'):g.plot(ax=ax,color=C['P']['COLORES_FALLAS'].get(name,'#555'),linestyle=':',linewidth=1.1,zorder=5)
 for rid,p in C['profiles'].items():
  geom=LineString(p[['x_utm','y_utm']].to_numpy());inside=geom.intersection(poly.geometry.union_all())
  if inside.is_empty:continue
  co='#0e4c82' if rid not in [4,5,6,7,8,9] else '#982264'
  for line in (inside.geoms if inside.geom_type=='MultiLineString' else [inside]):
   if line.geom_type=='LineString':xy=np.array(line.coords)/1000;ax.plot(xy[:,0],xy[:,1],c=co,lw=1.05,zorder=6)
  if C['P']['ROTULAR_RIOS_MAPA']:
   middle=inside.interpolate(.5,normalized=True);ax.annotate(f'R{rid:02d}',(middle.x/1000,middle.y/1000),fontsize=7,color=co,xytext=(3,0),textcoords='offset points')
 x0=b[0]/1000;y0=b[1]/1000;bar=C['P']['BARRA_ESCALA_KM'];xr=(b[2]-b[0])/1000;yr=(b[3]-b[1])/1000
 bx=x0+.06*xr;by=y0+.04*yr;ax.plot([bx,bx+bar],[by,by],c='black',lw=2);ax.text(bx+bar/2,by+.015*yr,f'{bar:g} km',ha='center',fontsize=9)
 ax.annotate('N',xy=(.94,.96),xytext=(.94,.86),xycoords='axes fraction',ha='center',arrowprops={'arrowstyle':'->','lw':1.5})
 ax.set_xlim(b[0]/1000,b[2]/1000);ax.set_ylim(b[1]/1000,b[3]/1000);ax.set_aspect('equal');ax.set_xlabel('Este UTM 18S (km)');ax.set_ylabel('Norte UTM 18S (km)');ax.grid(alpha=.15)
 fig.suptitle(f'Red de drenaje del polígono de estudio\nÁrea aportante ≥ {amin:g} km² · {len(edges):,} segmentos'.replace(',','.'))
 ax.legend(handles=[Line2D([],[],c=C['P']['COLOR_RED'],lw=1,label=f'Red A ≥ {amin:g} km²'),Line2D([],[],c='#982264',lw=1.5,label='Cauces con calibración K y τ'),Line2D([],[],c='#0e4c82',lw=1,label='Otros cauces principales'),Line2D([],[],c='#555',ls=':',label='Fallas CHAF')],loc='lower right',fontsize=8)
 guardar(C,fig,cat,'Mapa_regional_red_05km2','Mapa regional de la red de drenaje',
    f'Polígono suministrado PoligonoFallas y red D8 recortada a su geometría. Se muestran todos los segmentos con área aportante ≥ {amin:g} km², sin filtrar por orden ni por cruce de falla. Fondo: relieve del DEM original. Los aportes se conservaron desde el cálculo hidrológico global y no se recalcularon sobre un DEM recortado.')
 guardar_json(C['datos']/'Control_mapa_regional.json',{'fuente':str(src),'capa':layer,'area_minima_km2':amin,'n_segmentos':len(edges),'min_area_mostrada_km2':float(edges.area_km2.min()),'poligono_area_km2':float(poly.area.sum()/1e6),'sin_submuestrear_geometria':True})

def reloca(C):
 cat='11_Reloca_cuatro_cortes';iniciar(C,cat)
 import reloca_cuatro_cortes_figuras as m
 root=C['datos']/'Reloca_cuatro_cortes';F=root/'Fuentes';F.mkdir(parents=True,exist_ok=True);(root/'Resultados').mkdir(exist_ok=True)
 files=[(C['run']/'K_oficiales_tres_PT.csv','K_oficiales_tres_PT.csv'),(C['run']/'Comparacion/K_tres_PT.csv','K_comparacion_tres_PT.csv'),
  (C['run']/'perfil_R06_utilizado.csv','perfil_R06_utilizado.csv'),(C['run']/'Trazabilidad.json','Trazabilidad_K.json'),
  (Path(C['P']['CSV_KP_ORIGINAL']),'kp_manuales_28_09.csv')]
 path=Path(C['P']['CAPA_CORTES_RELOCA']);files.extend((p,p.name) for p in path.parent.glob(path.stem+'.*') if p.is_file())
 for src,name in files:shutil.copy2(src,F/name)
 guardar_json(F/'Manifest_fuentes.json',{'corrida_K':str(C['run']),'capa_cortes':str(path),'fuentes_sha256':{name:sha(F/name) for src,name in files}})
 m.ESTILO_FIGURAS=C['P'];summary,controls=m.calcular(root)
 summary.to_csv(C['out']/cat/'Resumen_cuatro_cortes.csv',index=False,encoding='utf-8-sig')
 for cid in [1,2,3,4]:
  q=summary[summary.id_corte.eq(cid)].set_index('MIS_donante');ref=q.loc['5e']
  detalle=(f'Tiempos con interceptos intermedios: B–corte {ref.tau_hasta_corte_intermedio_ka:.1f} ka; '
   f'total al KP con donante 5e {q.loc["5e"].tau_KP_intermedios_ka:.1f} ka y con donante 5c {q.loc["5c"].tau_KP_intermedios_ka:.1f} ka. '
   f'Conserva el ancla de {int(ref.n_PT_Reloca_preservados)} de los 3 Pₜ. ')
  for kind in ['chi_altura','tau_altura']:
   name=f'Corte_{cid}_{kind}'
   for ext in C['P']['FORMATOS']:
    source=root/'Resultados'/(name+'.'+ext)
    if source.exists():shutil.copy2(source,C['out']/cat/source.name)
   f=C['out']/'Catalogo_figuras.json';entries=json.loads(f.read_text(encoding='utf-8'))
   entries.append({'categoria':cat,'nombre':name,'titulo':f'Reloca corte {cid} '+('chi altura' if kind=='chi_altura' else 'tau altura'),'rid':6,'mis':None,
    'png':cat+'/'+name+'.png','formatos':C['P']['FORMATOS'],
    'leyenda':detalle+f'Corte {cid} como cambio alternativo de K: Reloca desde N04 hasta el corte y Curanilahue hasta el KP. Colores separados por MIS; líneas para los Pₜ a edad nominal, sin rellenos. Los cuatro cortes son alternativas independientes. La conservación de τ(Pₜ)=T no valida por sí misma el K transferido.',
    'parametros_sha256':hashlib.sha256(json.dumps(C['P'],sort_keys=True).encode()).hexdigest()})
   guardar_json(f,entries)
 return summary

def mostrar(C,cat):
 from IPython.display import display,HTML
 entries=json.loads((C['out']/'Catalogo_figuras.json').read_text(encoding='utf-8'));entries=[e for e in entries if e['categoria']==cat]
 print(f'{CATEGORIAS[cat]}: {len(entries)} figuras guardadas en {C["out"]/cat}')
 if C['P']['MOSTRAR_EN_NOTEBOOK']:
  for e in entries:
   if e['rid'] is None or e['rid'] in C['P']['RIOS_VISTA']:
    from IPython.display import Image
    display(Image(filename=str(C['out']/e['png']),width=C['P']['ANCHO_VISTA_PX']))

def informe(C):
 entries=json.loads((C['out']/'Catalogo_figuras.json').read_text(encoding='utf-8'))
 entries.sort(key=lambda e:(e['categoria'],e['rid'] is None,e['rid'] or 0,e['nombre']))
 missing=[e['png'] for e in entries if not (C['out']/e['png']).exists()]
 assert not missing,missing
 counts=pd.Series([e['categoria'] for e in entries]).value_counts().to_dict()
 guardar_json(C['out']/'Catalogo_figuras.json',entries)
 rows=['<!doctype html><html lang="es"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Informe de figuras fluviales</title>',
 '<style>body{font:16px system-ui;max-width:1500px;margin:25px auto;padding:20px;color:#183247}h1,h2{color:#111}nav{display:flex;flex-wrap:wrap;gap:12px;margin:25px 0}a{color:#1764a0}details{margin:20px 0;border-top:1px solid #ccc;padding:15px 0}summary{font-weight:600;font-size:19px;cursor:pointer}figure{margin:25px 0;padding:15px 0;border-bottom:1px solid #ddd}table{border-collapse:collapse;font-size:13px}td,th{padding:7px;border:1px solid #ccd3da}th{background:#e9f0f5}img{display:block;max-width:100%;max-height:950px;margin:auto}figcaption{font-size:14px;line-height:1.5;max-width:1100px;margin:15px auto}button,select{padding:7px;font:inherit}.controls{position:sticky;top:0;background:white;padding:12px 0;border-bottom:1px solid #ddd;z-index:3}@media print{details{display:block}summary{display:none}.controls,nav{display:none}figure{break-inside:avoid}}</style>',
 '<h1>Informe de figuras de geomorfología fluvial</h1><p>20 de octubre de 2026 · θ = 0,40 · A₀ = 10⁶ m² · n = 1</p>',
 f'<p>Se reúnen {len(entries)} figuras: perfiles y ksn de 18 cauces, calibraciones de K de seis cauces, perfiles completos hasta la costa y τ de los ríos con K propio o transferido, cuatro cortes alternativos de Reloca y la red del polígono regional desde {C["P"]["AREA_MINIMA_RED_KM2"]:g} km² de área aportante. Los ocho casos oficiales B/MIS conservan tres etiquetas Pₜ, equivalentes a 18 posiciones físicas contadas por caso.</p>',
 '<p><b>Cronología:</b> edades asignadas de highstand en Tabla 2 de Jara-Muñoz et al. (2015): MIS 5e = 125 ± 5 ka y MIS 5c = 106 ± 5 ka. La geometría de terrazas proviene del catálogo recibido del autor. Los promedios norte–sur se calculan dentro de cada MIS; los segmentos unen directamente los promedios de cada río, sin regresión. U de George es aparente; U del método de pendientes se impone desde la terraza.</p>',
 '<p><b>Alcance:</b> las fallas punteadas corresponden a intersecciones con trazas CHAF. La cobertura geológica es incompleta y las franjas indican la cartografía disponible, sin verificar la roca del lecho. Los KP proceden del CSV manual y se proyectan por coordenadas al perfil actual. Chanco conserva su figura τ, con aviso de que no hay KP manual registrado. Los tiempos son escenarios de propagación, no edades tectónicas directas.</p>',
 '<p><a href="Informe_figuras_20261020.pdf">Abrir informe PDF</a> · <a href="Parametros_figuras.json">Parámetros usados</a></p>',
 '<div class="controls"><label>Río visible <select id="rio"><option value="todos">Todos</option>'+''.join(f'<option value="{r}">R{r:02d} · {escape(C["names"][r])}</option>' for r in sorted(C['names']))+'</select></label> <button onclick="document.querySelectorAll(\'details\').forEach(x=>x.open=true)">Abrir secciones</button> <button onclick="document.querySelectorAll(\'details\').forEach(x=>x.open=false)">Cerrar secciones</button></div>',
 '<nav>'+''.join(f'<a href="#{cat}">{cat.split("_")[0]} · {escape(title)}</a>' for cat,title in CATEGORIAS.items())+'</nav>']
 for cat,title in CATEGORIAS.items():
  rows.append(f'<details id="{cat}"><summary>{escape(title)} · {counts.get(cat,0)} figuras</summary>')
  if cat=='13_Tablas_K_tau':
   alturas=pd.read_csv(C['out']/'Tabla_KP_altura_edades_MIS5e.csv',sep=';',decimal=',')
   rows.append('<h3>Altura de KP y tiempos con K MIS 5e</h3>'+alturas.to_html(index=False,na_rep='—',float_format=lambda x:f'{x:.2f}'))
   rows.append('<p><a href="Tabla_KP_altura_edades_MIS5e.csv">Descargar alturas y edades MIS 5e</a></p>')
   t=pd.read_csv(C['datos']/'K_tau_atributos/Resumen_K_promedio_tau.csv')
   cols=['rid','rio','mis','fuente_K','origen_tau','K_George','K_pendiente','tau_George_ka','tau_pendiente_ka','estado']
   rows.append('<p>K en m⁰·²/año; τ en ka. Medias separadas por MIS y método. τ se calcula con el K promedio. R01 queda sin cálculos. Las transferencias usan el extremo inferior del perfil receptor como origen. Reloca: aquí K uniforme; los cuatro cortes alternativos están en la sección 11.</p>')
   rows.append('<div style="overflow-x:auto">'+t[cols].to_html(index=False,na_rep='—',float_format=lambda x:f'{x:.6g}')+'</div>')
   rows.append('<p><a href="../Fuentes_vigentes_figuras_20261020/K_tau_atributos/Resumen_K_promedio_tau.csv">Descargar tabla CSV</a> · <a href="../Fuentes_vigentes_figuras_20261020/K_tau_atributos/kp_manuales_K_tau.gpkg">Capa de KP y escenarios GPKG</a></p>')
  if cat=='11_Reloca_cuatro_cortes':
   t=pd.read_csv(C['out']/cat/'Resumen_cuatro_cortes.csv')
   t=t[['id_corte','MIS_donante','tau_hasta_corte_intermedio_ka','tau_corte_KP_intermedio_ka','tau_KP_intermedios_ka','n_PT_Reloca_preservados']]
   t.columns=['Corte','MIS donante','τ B–corte (ka)','τ corte–KP (ka)','τ total KP (ka)','Pₜ conservados de 3']
   rows.append('<div style="overflow-x:auto">'+t.to_html(index=False,float_format=lambda x:f'{x:.2f}')+'</div>')
  for e in entries:
   if e['categoria']!=cat:continue
   stem=str(Path(e['png']).with_suffix('')).replace('\\','/')
   links=' · '.join(f'<a href="{escape(stem+"."+ext)}">{ext.upper()}</a>' for ext in e['formatos'] if (C['out']/(stem+'.'+ext)).exists())
   rows.append(f'<figure data-rid="{e["rid"] if e["rid"] is not None else "global"}"><h3>{escape(e["titulo"])}</h3><a href="{escape(e["png"])}"><img loading="lazy" src="{escape(e["png"])}" alt="{escape(e["titulo"])}"></a><figcaption>{escape(e["leyenda"])}<br>{links}</figcaption></figure>')
  rows.append('</details>')
 rows+=['<p>Referencias metodológicas: Perron y Royden (2013), doi:10.1002/esp.3302; Royden y Perron (2013), doi:10.1002/jgrf.20031; Jara-Muñoz et al. (2015), doi:10.1016/j.quascirev.2015.01.005.</p>',
 '<script>const filter=document.getElementById("rio");filter.onchange=()=>{document.querySelectorAll("figure").forEach(f=>f.hidden=filter.value!=="todos"&&f.dataset.rid!==filter.value&&f.dataset.rid!=="global")};function showHash(){let id=decodeURIComponent(location.hash.slice(1));let el=document.getElementById(id);if(el&&el.tagName==="DETAILS")el.open=true}window.addEventListener("hashchange",showHash);showHash();</script></html>']
 (C['out']/'Informe_figuras.html').write_text('\n'.join(rows),encoding='utf-8')
 if C['P']['GENERAR_PDF']:
  script=Path(__file__).with_name('informe_pdf_figuras_20261020.py')
  subprocess.run([C['P']['PYTHON_PDF'],str(script),str(C['out'])],check=True)
 report={'fecha_utc':datetime.now(timezone.utc).isoformat(),'n_figuras':len(entries),'por_categoria':counts,'archivos_verificados':True,
    'theta':.40,'n_rios_perfiles':18,'n_rios_calibrados':6,'PDF_generado':C['P']['GENERAR_PDF']}
 guardar_json(C['out']/'Verificacion_informe.json',report)
 print(json.dumps(report,ensure_ascii=False,indent=2),flush=True)
 return C['out']/'Informe_figuras.html'

_regional_drenaje=regional
def regional(C):
 _regional_drenaje(C)
 from revision_figuras_20261020 import regional_ksn
 regional_ksn(C)

def resumen(C):
 from revision_figuras_20261020 import resumen as generar
 return generar(C)

# Current display profiles and frozen calibration profiles have different roles.
CATEGORIAS={'00_Perfiles_costa_planicie':'Perfiles completos hasta la costa y recorte de planicie aluvial',**CATEGORIAS}
_cargar_calibracion_vigente=cargar
def cargar(P):
 from figuras_vigentes_20261020 import aplicar_perfiles
 return aplicar_perfiles(_cargar_calibracion_vigente(P))
_metodo_pendiente_calibracion=metodo_pendiente
def metodo_pendiente(C):
 from figuras_vigentes_20261020 import contexto_calibracion
 return _metodo_pendiente_calibracion(contexto_calibracion(C))
_metodo_george_calibracion=metodo_george
def metodo_george(C):
 from figuras_vigentes_20261020 import contexto_calibracion
 return _metodo_george_calibracion(contexto_calibracion(C))
_residuos_calibracion=residuos_y_planta
def residuos_y_planta(C):
 from figuras_vigentes_20261020 import contexto_calibracion
 return _residuos_calibracion(contexto_calibracion(C))
def perfiles_costa(C):
 from figuras_vigentes_20261020 import perfiles_costa as generar
 return generar(C)
def tau(C):
 from figuras_vigentes_20261020 import tau as generar
 return generar(C)
def tabla_alturas(C):
 from figuras_vigentes_20261020 import tabla_alturas as generar
 return generar(C)
