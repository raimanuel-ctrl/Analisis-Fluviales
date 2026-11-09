"""Revisión de tesis: cartografía ksn y atributos K/tau trazables."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
from datetime import datetime,timezone
import csv,json,shutil,hashlib
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from shapely.geometry import Point

def regional_ksn(C):
 import red_fallas as rf
 import figuras_memoria_20261020 as f
 cat='12_Mapa_regional_red';src=Path(C['P']['RED_REGIONAL'])
 poly=gpd.read_file(src,layer='poligono').to_crs(32718);geom=poly.geometry.union_all()
 edges=gpd.read_file(src,layer=C['P']['CAPA_RED'],where=f"area_km2 >= {C['P']['AREA_MINIMA_RED_KM2']}").to_crs(32718)
 graph=C['base']/'Red_drenaje_poligono_20261019/cache/grafo_06cbe8fc3f23a21dd9.npz'
 with np.load(graph) as z:g={k:z[k] for k in z.files}
 assert float(g['theta'])==.4 and float(g['A0'])==1e6
 ids=edges.node_u.to_numpy(int);assert np.array_equal(edges.gid_u,g['gid'][ids])
 cache=C['datos']/'Red_ksn_500_1000';cache.mkdir(exist_ok=True)
 controls=[];key=hashlib.sha256((f.sha(graph)+geom.wkb_hex).encode()).hexdigest()
 for W in C['P']['VENTANAS_KSN_REGIONAL_M']:
  file=cache/f'ksn_{int(W)}m.npz';meta=file.with_suffix('.json')
  expected={'grafo_poligono_sha256':key,'ventana_m':W,'min_n':C['P']['MIN_NODOS_KSN'],'min_fraction':C['P']['FRACCION_VENTANA_MIN'],'area_minima_visible_km2':float(C['P']['AREA_MINIMA_RED_KM2']),'gid_u_sha256':hashlib.sha256(edges.gid_u.to_numpy().tobytes()).hexdigest()}
  valid=file.exists() and meta.exists() and all(json.loads(meta.read_text()).get(k)==v for k,v in expected.items())
  if valid:
   print(f'ksn regional {W} m desde caché verificada.',flush=True)
   with np.load(file) as z:values={k:z[k] for k in z.files}
  else:
   print(f'Calculando ksn regional independiente: ventana {W} m.',flush=True)
   fields,inside,info=rf.calcular_ksn(g,geom,window=W,min_n=C['P']['MIN_NODOS_KSN'],min_fraction=C['P']['FRACCION_VENTANA_MIN'])
   values={k:a[ids] for k,a in fields.items()};np.savez_compressed(file,**values)
   f.guardar_json(meta,{**expected,**info,'n_segmentos':len(edges),'gid_u_sha256':hashlib.sha256(edges.gid_u.to_numpy().tobytes()).hexdigest()})
  edges[f'ksn{int(W)}']=values['ksn'];edges[f'R2_{int(W)}']=values['R2']
  controls.append({**expected,'n_segmentos':len(edges),'n_validos':int(np.isfinite(values['ksn']).sum()),'n_invalidos':int((~np.isfinite(values['ksn'])).sum())})
 vals=np.concatenate([edges[f'ksn{int(w)}'].dropna().to_numpy() for w in C['P']['VENTANAS_KSN_REGIONAL_M']])
 lim=C['P'].get('KSN_REGIONAL_LIM') or [0,float(np.quantile(vals,.98))]
 norm=Normalize(*lim,clip=True);bounds=poly.total_bounds;pad=max(bounds[2:]-bounds[:2])*.015
 b=tuple(bounds+np.array([-pad,-pad,pad,pad]));_,hs,bb=f.relieve(C,b,maxside=1600)
 segments=[];indexes=[]
 for i,line in enumerate(edges.geometry):
  for part in (line.geoms if line.geom_type=='MultiLineString' else [line]):
   segments.append(np.array(part.coords)/1000);indexes.append(i)
 indexes=np.array(indexes)
 for W in C['P']['VENTANAS_KSN_REGIONAL_M']:
  fig,ax=plt.subplots(figsize=(11.5,10),layout='constrained')
  if C['P']['MOSTRAR_RELIEVE_MAPA']:ax.imshow(hs,extent=np.array(bb)[[0,2,1,3]]/1000,origin='upper',cmap='gray',alpha=.45,vmin=0,vmax=1)
  ax.add_collection(LineCollection(segments,colors='#b6bdc4',linewidths=.38,rasterized=True,zorder=2))
  v=edges[f'ksn{int(W)}'].to_numpy()[indexes];ok=np.isfinite(v)
  lc=LineCollection([segments[i] for i in np.flatnonzero(ok)],cmap=C['P']['KSN_REGIONAL_CMAP'],norm=norm,linewidths=.7,rasterized=True,zorder=3)
  lc.set_array(v[ok]);ax.add_collection(lc)
  pk=poly.copy();pk.geometry=pk.geometry.scale(xfact=.001,yfact=.001,origin=(0,0));pk.boundary.plot(ax=ax,color='#222',lw=1)
  if C['P']['MOSTRAR_FALLAS_MAPA']:
   fk=C['faults'].copy();fk.geometry=fk.geometry.scale(xfact=.001,yfact=.001,origin=(0,0));fk.plot(ax=ax,color='#222',ls=':',lw=1,zorder=5)
  for rid,p in C['profiles'].items():
   mid=Point(float(p.x_utm.iloc[len(p)//2]),float(p.y_utm.iloc[len(p)//2]))
   if geom.covers(mid):ax.annotate(f'R{rid:02d}',(mid.x/1000,mid.y/1000),fontsize=8,bbox={'fc':'white','ec':'none','alpha':.8},zorder=7)
  cbar=fig.colorbar(lc,ax=ax,shrink=.65,pad=.02,extend='max');cbar.set_label('ksn (m⁰·⁸) · escala común a ambas ventanas')
  ax.set(xlim=(b[0]/1000,b[2]/1000),ylim=(b[1]/1000,b[3]/1000),xlabel='Este UTM 18S (km)',ylabel='Norte UTM 18S (km)');ax.set_aspect('equal');ax.grid(alpha=.12)
  bx=b[0]/1000+.05*(b[2]-b[0])/1000;by=b[1]/1000+.04*(b[3]-b[1])/1000
  bar=C['P']['BARRA_ESCALA_KM'];ax.plot([bx,bx+bar],[by,by],c='k',lw=2);ax.text(bx+bar/2,by+1,f'{bar:g} km',ha='center')
  ax.annotate('N',xy=(.95,.95),xytext=(.95,.85),xycoords='axes fraction',ha='center',arrowprops={'arrowstyle':'->'})
  from matplotlib.lines import Line2D
  ax.legend(handles=[Line2D([],[],color='#222',ls=':',label='Fallas CHAF'),Line2D([],[],color='#b6bdc4',label='Sin ajuste válido')],loc='upper left',fontsize=8,framealpha=.9)
  fig.suptitle(f'ksn regional · ventana {W:g} m · θ = 0,40\nÁrea aportante visible ≥ {C["P"]["AREA_MINIMA_RED_KM2"]:g} km²')
  f.guardar(C,fig,cat,f'Mapa_regional_ksn_{int(W)}m',f'ksn regional · ventana {W:g} m',
   f'OLS altura–χ con intercepto libre en ventana de {W:g} m, calculada sobre el grafo D8 completo con contexto desde 0,1 km²; se muestran todos los segmentos ≥0,5 km² del polígono. θ=0,40 y A₀=10⁶ m². Gris: ajuste inválido o ventana insuficiente; no se reemplaza por cero. Las dos ventanas comparten escala; saturación al percentil 98 combinado, conservando los valores completos en la capa. Fallas CHAF punteadas.')
 dest=cache/'Red_ksn_500_1000.gpkg';edges.to_file(dest,layer='red_ksn_500_1000',driver='GPKG')
 f.guardar_json(cache/'Control_ventanas.json',{'ventanas':controls,'limites_color':lim,'sin_submuestreo':True,'grafo':str(graph),'capa':str(dest)})
 print('Dos ventanas regionales guardadas.',flush=True)

def resumen(C):
 import figuras_memoria_20261020 as f
 cat='13_Tablas_K_tau';f.iniciar(C,cat)
 out=C['datos']/'K_tau_atributos';out.mkdir(exist_ok=True)
 rows=[];long=[];assign={};donors={int(k):int(v) for k,v in C['P']['DONANTES_EXPERIMENTALES'].items()}
 for rid in range(1,19):
  p=C['profiles'][rid];kp=C['kp'][C['kp'].rid.eq(rid)];own=C['k'][C['k'].rid.eq(rid)]
  source=rid if len(own) else donors.get(rid)
  if rid==1:source=None
  if source is None:
   rows.append({'rid':rid,'rio':C['names'][rid],'mis':'—','fuente_K':'Sin cálculo','origen_tau':'—','K_George':np.nan,'K_pendiente':np.nan,'tau_George_ka':np.nan,'tau_pendiente_ka':np.nan,'estado':'Excluido por instrucción'});continue
  for mis,ng in C['ns'][C['ns'].rid.eq(source)].groupby('mis'):
   mean=ng.iloc[0];base=C['k'][C['k'].rid.eq(source)&C['k'].mis.eq(mis)].iloc[0]
   chi0=float(np.interp(base.s_B_m,p.distancia_m,p.chi_m)) if source==rid else float(p.chi_m.iloc[0])
   s0=float(base.s_B_m) if source==rid else float(p.distancia_m.iloc[0])
   r={'rid':rid,'rio':C['names'][rid],'mis':mis,'fuente_K':f'R{source:02d}','rid_fuente_K':source,'experimental':source!=rid,'origen_tau':f'B {base.id_B_corto}' if source==rid else 'Extremo inferior del perfil receptor','chi_origen_m':chi0,'s_origen_m':s0,'theta':.4,'n_PT_distintos':int(mean.n_PT_distintos),'K_George':float(mean.K_George),'K_pendiente':float(mean.K_pendiente),'tau_George_ka':np.nan,'tau_pendiente_ka':np.nan,'estado':'Sin KP manual' if kp.empty else ('Transferencia experimental' if source!=rid else 'Calibración propia')}
   for point in kp.itertuples():
    assert point.chi_m>=chi0,('KP aguas abajo del origen',rid,mis)
    attrs=assign.setdefault(point.id_kp,{'rid':rid,'k_src':source,'exp_K':int(source!=rid),'theta_tau':.4,'chi040_m':point.chi_m})
    attrs[f'chi0_{mis}']=chi0;attrs[f's0_{mis}_m']=s0
    for method,prefix in [('George','G'),('pendiente','P')]:
     kval=r['K_'+method];tau=(point.chi_m-chi0)/(kval*1e6**.4*1000)
     assert np.isfinite(tau) and tau>=0
     r['tau_'+method+'_ka']=float(tau)
     attrs[f'K{prefix}_{mis}']=kval;attrs[f't{prefix}_{mis}_ka']=tau
     long.append({'rid':rid,'rio':C['names'][rid],'id_kp':point.id_kp,'mis':mis,'metodo':method,'K_promedio':kval,'tau_KP_ka':tau,'rid_fuente_K':source,'experimental':source!=rid,'origen_tau':r['origen_tau'],'chi_origen_m':chi0,'chi_KP_m':point.chi_m,'theta':.4,'A0_m2':1e6,'n_PT_distintos':int(mean.n_PT_distintos),'separacion_cauce_m':point.separacion_cauce_m})
   rows.append(r)
 table=pd.DataFrame(rows);times=pd.DataFrame(long)
 table.to_csv(out/'Resumen_K_promedio_tau.csv',index=False,encoding='utf-8-sig');times.to_csv(out/'Tau_por_K_promedio_escenarios.csv',index=False,encoding='utf-8-sig')
 # Tables are last in both the offline HTML and the PDF, and have editable originals.
 for kind,sel in [('Calibraciones_propias',table[table.rid.between(4,9)]),('Transferencias_experimentales',table[~table.rid.between(4,9)])]:
  values=[]
  for v in sel.itertuples():
   num=lambda x:'—' if pd.isna(x) else f'{x:.2f}'
   values.append([f'R{v.rid:02d}',v.rio.replace('Estero de ','').replace('Estero ','').replace('Rio ',''),v.mis,v.fuente_K,num(v.K_George*1e5),num(v.K_pendiente*1e5),num(v.tau_George_ka),num(v.tau_pendiente_ka)])
  fig,ax=plt.subplots(figsize=(15,max(5,1.8+.34*len(sel))));ax.axis('off')
  t=ax.table(cellText=values,colLabels=['ID','Río','MIS','Fuente K','K George\n×10⁻⁵','K pendiente\n×10⁻⁵','τ George\nka','τ pendiente\nka'],loc='center',cellLoc='center',colWidths=[.05,.24,.06,.08,.13,.14,.15,.15]);t.auto_set_font_size(False);t.set_fontsize(10);t.scale(1,1.65)
  for (i,j),cell in t.get_celld().items():
   cell.set_edgecolor('#d5dde3')
   if i==0:cell.set_facecolor('#203e53');cell.set_text_props(color='white',weight='bold');cell.set_height(cell.get_height()*1.5)
   elif i%2==0:cell.set_facecolor('#f0f4f7')
  ax.set_title(('K y τ · calibraciones propias' if kind.startswith('Calib') else 'K y τ · transferencias experimentales')+'\nK en m⁰·²/año; medias por MIS y método · θ=0,40',fontsize=14,pad=22)
  fig.tight_layout()
  f.guardar(C,fig,cat,kind,kind.replace('_',' '),'K promedio aritmético de Pₜ físicos distintos dentro de cada río, MIS y método. τ se calcula con ese K promedio, no promediando tiempos. Propios: origen B de cada MIS; transferidos: extremo inferior del perfil receptor, τ=0. R02/R03 reciben R04; R10–R18 reciben R09. R01 excluido y Chanco sin KP manual. Los tiempos de Reloca aquí usan K uniforme; los cuatro cambios de K se conservan como experimentos independientes en la sección 11.')
 mapping={'rid':'id del río receptor','k_src':'río fuente de K','exp_K':'1=transferencia experimental; 0=K propio','theta_tau':'theta del cálculo de tau','chi040_m':'chi KP recalculado con theta 0.40'}
 for mis in ['5e','5c']:
  mapping.update({f'chi0_{mis}':'chi del origen de tau; '+mis,f's0_{mis}_m':'distancia al origen en perfil receptor; '+mis})
  for m,prefix in [('George','G'),('pendiente','P')]:mapping.update({f'K{prefix}_{mis}':f'K promedio {m}, MIS {mis}, m^0.2/año',f't{prefix}_{mis}_ka':f'tau con K promedio {m}, MIS {mis}, ka'})
 target=Path(C['P']['CSV_KP_ORIGINAL']);raw=target.read_bytes()
 with target.open(encoding='utf-8-sig',newline='') as stream:
  reader=csv.DictReader(stream,delimiter=';');original_fields=reader.fieldnames;original=list(reader)
 valid_manual=[r for r in original if r.get('x_utm') and r.get('y_utm')]
 assert {r['id_kp'] for r in valid_manual}==set(C['kp'].id_kp),'Las selecciones del CSV original difieren de la corrida consolidada.'
 for r in valid_manual:
  kp=C['kp'].set_index('id_kp').loc[r['id_kp']]
  assert abs(float(r['x_utm'].replace(',','.'))-kp.x_original_m)<.1
  assert abs(float(r['y_utm'].replace(',','.'))-kp.y_original_m)<.1
 # Preserve all existing lexical values; append/replace derived fields only.
 for row in original:
  rid=int(float(row['id'].replace(',','.')));a=assign.get(row.get('id_kp'),{})
  for name in mapping:
   val=a.get(name);row[name]='' if val is None else (f'{val:.2f}' if name.startswith(('tG_','tP_')) and name.endswith('_ka') else str(val)).replace('.',',')
  if rid==1:assert all(row[n]=='' for n in mapping)
 fields=list(dict.fromkeys([*original_fields,*mapping]))
 enriched=out/'kp_manuales_con_K_tau.csv'
 with enriched.open('w',encoding='utf-8-sig',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=fields,delimiter=';');writer.writeheader();writer.writerows(original)
 f.guardar_json(out/'Diccionario_campos.json',{'campos':mapping,'promedio':'P_T físicos distintos; MIS separados','origen_calibrado':'B propio por MIS','origen_transferencia':'extremo inferior del perfil receptor, no B del donante','fuente_manual':str(target),'CSV_consolidado_congelado':'Se preserva sin cambios: sólo aporta las selecciones XY de esta corrida','R01':'Sin cálculos','R07':'Sin KP registrado','Reloca':'Tau uniforme aquí; cuatro cortes en sección 11'})
 # Export the same manual points, never the 31 unrelated old candidate locations.
 geo=[]
 for row in original:
  if not row.get('x_utm') or not row.get('y_utm'):continue
  parse=lambda s:float(s.replace(',','.'))
  rec={'rio_id':int(parse(row['id'])),'id_kp':row['id_kp'],'rio':row['rio']}
  for name in mapping:rec[name]=parse(row[name]) if row[name] else None
  rec['geometry']=Point(parse(row['x_utm']),parse(row['y_utm']));geo.append(rec)
 points=gpd.GeoDataFrame(geo,crs=32718)
 points.to_file(out/'kp_manuales_K_tau.gpkg',layer='kp_manuales',driver='GPKG')
 times_geo=times.merge(points[['id_kp','geometry']],on='id_kp',validate='many_to_one')
 gpd.GeoDataFrame(times_geo,geometry='geometry',crs=32718).to_file(out/'kp_manuales_K_tau.gpkg',layer='tau_escenarios',driver='GPKG')
 shp=target.with_suffix('.shp')
 if C['P']['ACTUALIZAR_KP_FUENTES']:
  stamp=datetime.now().strftime('%Y%m%d_%H%M%S');backup=out/'Respaldos'/stamp;backup.mkdir(parents=True)
  shutil.copy2(target,backup/target.name)
  for file in target.parent.glob(target.stem+'.*'):
   if file.suffix.lower() in ['.shp','.dbf','.shx','.prj','.cpg','.qmd','.qix']:shutil.copy2(file,backup/file.name)
  shutil.copy2(enriched,target);points.to_file(shp,driver='ESRI Shapefile',encoding='UTF-8')
  with target.open(encoding='utf-8-sig',newline='') as stream:after=list(csv.DictReader(stream,delimiter=';'))
  before=list(csv.DictReader(raw.decode('utf-8-sig').splitlines(),delimiter=';'))
  preserved=[k for k in original_fields if k not in mapping]
  assert len(before)==len(after)==18
  assert all(all(a[k]==b[k] for k in preserved) for a,b in zip(before,after))
  check=gpd.read_file(shp);assert len(check)==17 and check.crs.to_epsg()==32718
  assert check.loc[check.rio_id.eq(1),list(mapping)].drop(columns=['rid'],errors='ignore').isna().all().all()
  f.guardar_json(out/'Control_actualizacion_fuentes.json',{'fecha_utc':datetime.now(timezone.utc).isoformat(),'csv_actualizado':str(target),'shp_manual':str(shp),'backup':str(backup),'n_filas_csv':18,'n_puntos_shp':17,'n_escenarios_tau':len(times),'columnas_originales_preservadas':preserved,'verificado_lectura_posterior':True,'sin_K_R01':True,'sin_KP_inventado_R07':True,'CSV_original_sha256_antes':hashlib.sha256(raw).hexdigest(),'CSV_sha256_despues':f.sha(target)})
 print(table[['rid','mis','K_George','K_pendiente','tau_George_ka','tau_pendiente_ka']].to_string(index=False),flush=True)
 print('Atributos CSV y SHP manual:',target,shp,flush=True)
 return table
