"""Perfiles y tiempos vigentes; conserva las calibraciones originales de K."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd,geopandas as gpd,rasterio
from scipy.integrate import cumulative_trapezoid
from shapely.geometry import Point
import matplotlib.pyplot as plt

def aplicar_perfiles(C):
 import figuras_memoria_20261020 as f
 import red_fallas as rf
 from base_dos_K import project
 C['profiles_calibracion']={r:p.copy() for r,p in C['profiles'].items()}
 C['kp_calibracion']=C['kp'].copy();C['crossings_calibracion']=C['crossings'].copy()
 manifest=json.loads(Path(C['P']['MANIFIESTO_PERFILES_VIGENTES']).read_text(encoding='utf-8'))
 C['manifest_vigente']=manifest
 overrides={};checks=[]
 for rid in [9,13,14,15]:
  src=manifest['perfiles'][str(rid)];p=pd.read_csv(src['archivo'],float_precision='round_trip');p['chi_m']=p[src['columna_chi']];p['area_m2']=p[src['columna_area']]
  calc=1e6**.4*cumulative_trapezoid(p.area_m2.to_numpy()**-.4,p.distancia_m,initial=0)
  err=float(np.max(np.abs(calc-(p.chi_m-p.chi_m.iloc[0]))));assert err<1e-6,(rid,err)
  p=p.drop(columns=['codigo_lito','descripcion_lito','color_lito','estado_lito','fuentes_lito'],errors='ignore')
  overrides[rid]=p;checks.append({'rid':rid,'perfil':src['archivo'],'sha256':f.sha(src['archivo']),'error_chi_m':err})
 lito_cols=['codigo_lito','descripcion_lito','color_lito','estado_lito','fuentes_lito'];pending={};missing={}
 for rid,p in overrides.items():
  reference=C['profiles_calibracion'][rid][['x_utm','y_utm',*lito_cols]].drop_duplicates(['x_utm','y_utm'])
  merged=p.merge(reference,on=['x_utm','y_utm'],how='left',sort=False,validate='many_to_one');assert len(merged)==len(p)
  overrides[rid]=merged;missing[rid]=merged.codigo_lito.isna()
  if missing[rid].any():pending[rid]=merged.loc[missing[rid]].drop(columns=lito_cols).reset_index(drop=True)
 assigned,_=rf.asignar_litologia(pending,C['base']/C['cfg']['geologia_integrada']) if pending else ({},None)
 for rid,p in assigned.items():overrides[rid].loc[missing[rid],lito_cols]=p[lito_cols].to_numpy()
 for rid,p in overrides.items():
  p['color_lito']=[C['P']['COLORES_LITO'].get(str(c),v if f.is_color_like(v) else '#8b949e') for c,v in zip(p.codigo_lito,p.color_lito)]
  stats=rf.ventana_ols(p.distancia_m,p.chi_m,p.elevacion_m,W=C['P']['VENTANA_KSN_M'],min_n=C['P']['MIN_NODOS_KSN'],min_fraction=C['P']['FRACCION_VENTANA_MIN'],theta=.4,A0=1e6)
  p['ksn_040_figuras']=stats['ksn'];p['R2_ksn_figuras']=stats['R2'];p['ventana_efectiva_m']=stats['span'];C['profiles'][int(rid)]=p
 manual=pd.read_csv(C['P']['CSV_KP_ORIGINAL'],sep=';',decimal=',',float_precision='round_trip');kps=[]
 C['names'].update({int(r.id):r.rio for r in manual.itertuples()})
 for r in manual.itertuples():
  if not np.isfinite([r.id,r.x_utm,r.y_utm]).all():continue
  rid=int(r.id);at=project(C['profiles'][rid],[r.x_utm,r.y_utm]);assert at['separacion_cauce_m']<.1
  assert abs(at['z_m']-r.elevacion_m)<.1
  kps.append({'rid':rid,'id_kp':r.id_kp,'theta_seleccion_original':r.theta_seleccion,**at})
 C['kp']=pd.DataFrame(kps);C['crossings']=rf.cruces_en_perfiles(C['profiles'],C['faults'])
 if C['crossings'].empty:C['crossings']=pd.DataFrame(columns=['perfil_id','f_id','falla','nombre_base','chi_cruce_m','dist_cruce_m'])
 else:
  C['crossings']['qx']=(C['crossings'].x_cruce*100).round();C['crossings']['qy']=(C['crossings'].y_cruce*100).round()
  C['crossings']=C['crossings'].drop_duplicates(['perfil_id','nombre_base','qx','qy']).drop(columns=['qx','qy'])
 for i,v in C['k'].iterrows():
  p=C['profiles'][int(v.rid)];chi0=float(np.interp(v.s_B_m,p.distancia_m,p.chi_m));C['k'].loc[i,'chi_B_tau_m']=chi0
  q=C['kp'][C['kp'].rid.eq(int(v.rid))]
  if len(q):
   for method in ['George','pendiente']:C['k'].loc[i,'tau_KP_'+method+'_ka']=(q.chi_m.iloc[0]-chi0)/(v['K_'+method]*1e6**.4*1000)
 for rid,p in C['profiles'].items():p.to_csv(C['datos']/f'R{rid:02d}_perfil_figuras.csv',index=False,encoding='utf-8-sig')
 C['kp'].to_csv(C['datos']/'KP_proyectados_theta040.csv',index=False,encoding='utf-8-sig');C['crossings'].to_csv(C['datos']/'Cruces_fallas_18_perfiles.csv',index=False,encoding='utf-8-sig');C['k'].to_csv(C['datos']/'Calibraciones_figuras.csv',index=False,encoding='utf-8-sig')
 f.guardar_json(C['datos']/'Control_perfiles_vigentes.json',{'perfiles_actualizados':checks,'K_recalibrado':False,'K_y_U_fuente':str(C['run']),'R14':'P14 como outlet y KP C2 conservado','manual_actual':C['P']['CSV_KP_ORIGINAL'],'n_KP':len(C['kp']),'theta':.4})
 print('Actualizados: R09 con cuenca P5 recortada; R13=P13; R14=P14 y C2; R15=primer P15. K conservados.',flush=True)
 return C

def contexto_calibracion(C):
 return {**C,'profiles':C['profiles_calibracion'],'kp':C['kp_calibracion'],'crossings':C['crossings_calibracion']}

def construir_costa(C):
 import figuras_memoria_20261020 as f
 O=C['datos']/'Costa';O.mkdir(exist_ok=True);B=C['base'];coastpath=Path(C['P']['CAPA_BOCAS_COSTA']);coast=gpd.read_file(coastpath).to_crs(32718).set_index('id')
 OFF={1:(0,1),2:(1,1),4:(1,0),8:(1,-1),16:(0,-1),32:(-1,-1),64:(-1,0),128:(-1,1)}
 controls=[]
 with rasterio.open(B/'Productos/Datos/v1/topografia/flow_direction_D8_12m_v1.tif') as flow,rasterio.open(C['dem']) as dem:
  fd=flow.read(1);H,W=fd.shape
  for rid in range(1,19):
   old=pd.read_csv(B/f'Productos/Datos/v1/perfiles/rio_{rid:02d}_v1.csv',float_precision='round_trip');r,c=flow.index(old.x_utm.iloc[0],old.y_utm.iloc[0]);seen=set();path=[]
   for k in range(100000):
    if not(0<=r<H and 0<=c<W) or (r,c) in seen:break
    seen.add((r,c));path.append((r,c));code=int(fd[r,c])
    if code not in OFF:break
    dr,dc=OFF[code];r+=dr;c+=dc
   rr,cc=np.array(path).T;xx,yy=rasterio.transform.xy(flow.transform,rr,cc);xy=np.c_[xx,yy];point=coast.loc[rid].geometry;dist=np.linalg.norm(xy-[point.x,point.y],axis=1);i=int(np.argmin(dist));i_boca=i
   current=C['profiles'][rid];curxy=current[['x_utm','y_utm']].iloc[0].to_numpy();dcur=np.linalg.norm(xy-curxy,axis=1);ic=int(np.argmin(dcur))
   # Rahue's adopted downstream extension passes its historical mouth marker.
   # Include that real D8 segment rather than projecting the adopted outlet sideways.
   if dcur[ic]<.1 and ic>i:i=ic
   lowerxy=xy[:i+1][::-1];s=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(lowerxy,axis=0),axis=1))];L=float(s[-1]);zz=np.array([v[0] for v in dem.sample(lowerxy)])
   assert float(dist[i_boca])<24,(rid,float(dist[i_boca]));assert np.isfinite(zz).all()
   lower=pd.DataFrame({'x_utm':lowerxy[:,0],'y_utm':lowerxy[:,1],'distancia_costa_m':s,'elevacion_m':zz,'recorte_original':True})
   upper=old[['x_utm','y_utm','elevacion_m']].copy();upper['distancia_costa_m']=old.distancia_m+L;upper['recorte_original']=False;full=pd.concat([lower.iloc[:-1],upper],ignore_index=True)
   dd=np.linalg.norm(full[['x_utm','y_utm']].to_numpy()-curxy,axis=1);j=int(np.argmin(dd));assert dd[j]<.1,(rid,float(dd[j]))
   row={'rid':rid,'boca_fuente':str(coastpath),'distancia_boca_traza_m':float(dist[i_boca]),'limite_recorte_original_desde_costa_m':L,'outlet_actual_desde_costa_m':float(full.distancia_costa_m.iloc[j]),'distancia_outlet_actual_traza_m':float(dd[j]),'referencia_costera_aguas_arriba_salida':i>i_boca,'n_nodos':len(full),'archivo':str(O/f'R{rid:02d}_perfil_costa.csv')}
   full.to_csv(row['archivo'],index=False,encoding='utf-8-sig');controls.append(row)
 pd.DataFrame(controls).to_csv(O/'Control_perfiles_costa.csv',index=False,encoding='utf-8-sig');return pd.DataFrame(controls)

def perfiles_costa(C):
 import figuras_memoria_20261020 as f
 cat='00_Perfiles_costa_planicie';f.iniciar(C,cat);control=construir_costa(C)
 for v in control.itertuples():
  rid=int(v.rid);p=pd.read_csv(v.archivo);current=C['profiles'][rid];s0=v.outlet_actual_desde_costa_m;L=v.limite_recorte_original_desde_costa_m
  fig,axes=plt.subplots(1,2,figsize=(14,6),layout='constrained')
  for ax in axes:
   if C['P']['MOSTRAR_PLANICIE'] and L>0:ax.axvspan(0,L/1000,fc=C['P']['COLOR_PLANICIE'],alpha=.16,label='Tramo inferior recortado originalmente (planicie)')
   ax.plot(p.distancia_costa_m/1000,f.suavizar(C,p.elevacion_m,p.distancia_costa_m),c='#87939a',lw=1.1,label='Perfil completo hasta referencia costera')
   ax.plot((current.distancia_m+s0)/1000,f.suavizar(C,current.elevacion_m,current.distancia_m),c=C['P']['COLOR_PERFIL'],lw=1.25,label='Perfil utilizado actualmente')
   if C['P']['MOSTRAR_LIMITE_RECORTE_ORIGINAL'] and L>0:ax.axvline(L/1000,c='#ae752d',ls='--',lw=1,label='Límite del recorte original')
   ax.axvline(s0/1000,c='#247947',ls=':',lw=1.3,label='Outlet adoptado')
   for kp in C['kp'][C['kp'].rid.eq(rid)].itertuples():
    dd=np.hypot(p.x_utm-kp.x_original_m,p.y_utm-kp.y_original_m);j=int(np.argmin(dd));assert dd.iloc[j]<.1
    ax.scatter(p.distancia_costa_m.iloc[j]/1000,kp.z_m,c=C['P']['COLOR_KP'],s=C['P']['TAMANO_KP'],marker=C['P']['SIMBOLOS']['KP'],zorder=8,ec='white',lw=.5)
    ax.annotate('KP',(p.distancia_costa_m.iloc[j]/1000,kp.z_m),xytext=(5,7),textcoords='offset points',fontsize=8)
   if C['P']['MOSTRAR_FALLAS']:
    for crossing in C['crossings'][C['crossings'].perfil_id.eq(rid)].itertuples():ax.axvline((crossing.dist_cruce_m+s0)/1000,c=C['P']['COLORES_FALLAS'].get(crossing.nombre_base,'#555'),ls=':',lw=.8,alpha=.7)
   ax.set(xlabel='Distancia aguas arriba desde referencia costera (km)',ylabel='Altura DEM (m)');ax.grid(alpha=.18)
  axes[0].set_xlim(0,float(p.distancia_costa_m.iloc[-1])/1000);axes[0].set_title('Perfil completo')
  xmax=min(float(p.distancia_costa_m.iloc[-1]),max(L,s0)+C['P']['MARGEN_DETALLE_COSTA_M']);q=p[p.distancia_costa_m.le(xmax)]
  axes[1].set_xlim(0,max(xmax/1000,.1));axes[1].set_ylim(float(q.elevacion_m.min())-5,float(q.elevacion_m.max())+max(12,float(q.elevacion_m.max()-q.elevacion_m.min())*.12));axes[1].set_title('Detalle del tramo inferior')
  axes[0].legend(fontsize=7.1,loc='upper left');fig.suptitle(f.titulo(C,rid)+' · perfil hasta la costa y recorte aluvial',fontsize=14)
  if L==0:axes[1].text(.03,.96,'El perfil original ya comienza en la referencia costera;\nno hay un recorte inferior adicional documentado.',transform=axes[1].transAxes,va='top',fontsize=8)
  caption='Perfil completo reconstruido siguiendo D8 desde el outlet original a la boca de la capa histórica, con alturas del DEM original. Sombreado naranja: tramo retirado al ubicar el outlet original aguas arriba; corresponde al recorte manual atribuido a planicie, no a una clasificación automática del sustrato. Azul: perfil vigente; línea verde: outlet adoptado; línea naranja: límite original. La extensión costera es descriptiva y no modifica el origen de τ ni los K.'
  if v.referencia_costera_aguas_arriba_salida:caption+=' En Rahue la salida vigente queda aguas abajo de la referencia histórica de boca; se conserva el tramo D8 real hasta la salida vigente.'
  f.guardar(C,fig,cat,f'R{rid:02d}_costa_recorte',f.titulo(C,rid)+' Costa y planicie aluvial',caption,rid)
 print('18 perfiles completos exportados; outlets vigentes coinciden con el trazado.',flush=True)

def tau(C):
 import figuras_memoria_20261020 as f
 cat='06_Tau_knickpoints';f.iniciar(C,cat);rows=[]
 for rid in range(2,19):
  p=C['profiles'][rid];own=C['k'][C['k'].rid.eq(rid)];source=rid if len(own) else int(C['P']['DONANTES_EXPERIMENTALES'][rid]);kp=C['kp'][C['kp'].rid.eq(rid)]
  fig,axes=plt.subplots(1,2,figsize=(13,6),sharey=True,layout='constrained')
  for ax,method in zip(axes,['George','pendiente']):
   for v in C['ns'][C['ns'].rid.eq(source)].itertuples():
    mis=str(v.mis);K=float(getattr(v,'K_'+method));base=C['k'][C['k'].rid.eq(source)&C['k'].mis.eq(mis)].iloc[0];s0=float(base.s_B_m) if source==rid else float(p.distancia_m.iloc[0]);chi0=float(np.interp(s0,p.distancia_m,p.chi_m))
    end=float(kp.s_m.max()) if len(kp) and C['P']['TAU_HASTA_KP'] else float(p.distancia_m.iloc[-1]);s=np.unique(np.r_[s0,p.loc[p.distancia_m.between(s0,end),'distancia_m'],end]);xx=np.interp(s,p.distancia_m,p.chi_m);zz=np.interp(s,p.distancia_m,p.elevacion_m)
    if source==rid and C['P']['MOSTRAR_ESCENARIOS_PT_TAU']:
     for j,c in enumerate(own[own.mis.eq(mis)].itertuples()):ax.plot((xx-chi0)/(getattr(c,'K_'+method)*1e6**.4*1000),zz,c=f.color(C,mis),lw=.7,alpha=.45,ls=[':', '--','-.'][j%3])
    ax.plot((xx-chi0)/(K*1e6**.4*1000),zz,c=f.color(C,mis),lw=1.4,label=f'MIS {mis} · K promedio R{source:02d}')
    for point in kp.itertuples():
     tk=(point.chi_m-chi0)/(K*1e6**.4*1000);ax.scatter(tk,point.z_m,c=f.color(C,mis),marker=C['P']['SIMBOLOS']['KP'],s=C['P']['TAMANO_KP'],ec='white',lw=.5,zorder=5)
     ax.annotate(f'{tk:.2f} ka',(tk,point.z_m),xytext=(-5,-20 if mis=='5e' else 10),textcoords='offset points',ha='right',fontsize=8,color=f.color(C,mis))
     rows.append({'rid':rid,'rio':C['names'][rid],'mis':mis,'metodo':method,'K_promedio':K,'id_kp':point.id_kp,'altura_KP_m':point.z_m,'tau_KP_ka':tk,'fuente_K':source,'chi_origen_m':chi0,'chi_KP_m':point.chi_m,'origen_tau':'B '+str(base.id_B_corto) if source==rid else C['manifest_vigente']['perfiles'][str(rid)]['nombre_origen']})
   ax.set(xlabel='τ desde B (ka)' if source==rid else 'τ desde outlet adoptado (ka)',title='K por '+method,xlim=(0,None));ax.grid(alpha=.2);ax.legend(fontsize=8,loc='lower right')
  axes[0].set_ylabel('Altura DEM (m)');fig.suptitle(f.titulo(C,rid)+' · altura–τ · θ = 0,40'+(' · K transferido de R'+str(source).zfill(2) if source!=rid else '')+('\nSin KP manual registrado' if kp.empty else ''))
  f.guardar(C,fig,cat,f'R{rid:02d}_altura_tau',f.titulo(C,rid)+' Tiempos vigentes de KP','Líneas continuas: K promedio de Pₜ físicos distintos, por MIS y método. Estrellas: KP manual vigente, etiquetado con τ en ka. Sin áreas rellenas. Las líneas finas opcionales conservan los escenarios individuales de Pₜ de las calibraciones. R09 usa su cuenca recortada P5 con K fijo; R13 parte de P13, R14 de P14 con KP C2, R15 del primer P15. Los tiempos transferidos usan el outlet propio como origen y conservan la identificación experimental.',rid)
 pd.DataFrame(rows).to_csv(C['datos']/'Tau_KP_dos_metodos.csv',index=False,encoding='utf-8-sig')

def tabla_alturas(C):
 import figuras_memoria_20261020 as f
 times=pd.read_csv(C['datos']/'K_tau_atributos/Resumen_K_promedio_tau.csv');records=[]
 for rid in range(1,19):
  q=C['kp'][C['kp'].rid.eq(rid)];t=times[times.rid.eq(rid)&times.mis.eq('5e')]
  records.append({'ID':f'R{rid:02d}','Rio':C['names'][rid],'Altura_KP_m':float(q.z_m.iloc[0]) if len(q) else np.nan,'Tau_George_5e_ka':float(t.tau_George_ka.iloc[0]) if len(t) else np.nan,'Tau_pendiente_5e_ka':float(t.tau_pendiente_ka.iloc[0]) if len(t) else np.nan})
 table=pd.DataFrame(records);table.to_csv(C['out']/'Tabla_KP_altura_edades_MIS5e.csv',sep=';',decimal=',',float_format='%.2f',index=False,encoding='utf-8-sig')
 (C['out']/'Tabla_KP_altura_edades_MIS5e.html').write_text('<meta charset="utf-8"><h2>Altura y τ de KP con K MIS 5e</h2>'+table.to_html(index=False,na_rep='—',float_format=lambda x:f'{x:.2f}'),encoding='utf-8')
 return table
