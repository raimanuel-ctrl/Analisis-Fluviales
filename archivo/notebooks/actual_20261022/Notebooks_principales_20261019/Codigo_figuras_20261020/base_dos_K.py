
import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,itertools
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

FACTOR=1e6**.4

def project(p,point):
    xy=p[['x_utm','y_utm']].to_numpy();segments=np.diff(xy,axis=0);lengths=np.linalg.norm(segments,axis=1)
    assert (lengths>0).all()
    s=np.r_[0,np.cumsum(lengths)];f=np.clip(((np.asarray(point)-xy[:-1])*segments).sum(axis=1)/lengths**2,0,1)
    snap=xy[:-1]+f[:,None]*segments;distance=np.linalg.norm(snap-point,axis=1);i=int(distance.argmin());where=s[i]+f[i]*lengths[i]
    if distance[i]>35:raise ValueError('Punto a más de 35 m del perfil: revisar geometría.')
    return {'x_original_m':float(point[0]),'y_original_m':float(point[1]),'x_perfil_m':float(snap[i,0]),'y_perfil_m':float(snap[i,1]),
            'separacion_cauce_m':float(distance[i]),'s_m':float(np.interp(where,s,p.distancia_m)),
            'chi_m':float(np.interp(where,s,p.chi_m)),'z_m':float(np.interp(where,s,p.elevacion_m))}

def escenario(r,d,cut,kp):
    cb,ct,ck=r.chi_B_m,cut['chi_m'],kp['chi_m'];kr,kd=float(r.K_m02_anio),float(d.K_m02_anio)
    if not cb<ct<ck:raise ValueError('Se requiere B < corte < KP a lo largo de χ.')
    tc=(ct-cb)/(kr*FACTOR*1000);dt=(ck-ct)/(kd*FACTOR*1000)
    pt=float(r.chi_PT_m)
    tpt=(min(pt,ct)-cb)/(kr*FACTOR*1000)+max(pt-ct,0)/(kd*FACTOR*1000)
    time_after=max(pt-ct,0)/(kd*FACTOR*1000)
    available=float(r.T_ka)-time_after
    possible=available>0
    return {'B_Reloca':r.id_B_corto,'PT_Reloca':r.escenario_PT,'id_Reloca':r.id_medicion,
        'MIS_Reloca':r.mis,'MIS_donante':d.mis,'B_donante':d.id_B_corto,'PT_donante':d.escenario_PT,'id_donante':d.id_medicion,
        'K_Reloca':kr,'K_donante':kd,'T_Reloca_ka':float(r.T_ka),'T_donante_ka':float(d.T_ka),
        'tau_corte_ka':tc,'tau_corte_KP_ka':dt,'tau_KP_ka':tc+dt,
        'tau_KP_un_K_Reloca_ka':(ck-cb)/(kr*FACTOR*1000),
        'PT_aguas_arriba_del_corte':bool(pt>ct),'tau_corte_PT_ka':time_after,'tau_PT_dos_K_ka':tpt,
        'residuo_tau_PT_menos_T_ka':tpt-float(r.T_ka),'tiempo_disponible_B_corte_ka':available,
        'K_donante_minimo_para_K_inicial_positivo':max(pt-ct,0)/(FACTOR*float(r.T_ka)*1000),
        'K_inicial_positivo_posible':possible,
        'K_inicial_ajustado_si_posible':((min(pt,ct)-cb)/(FACTOR*available*1000) if possible else np.nan),
        'estatus':'escenario de transferencia; no compatible con calibración original' if abs(tpt-float(r.T_ka))>1e-8 else 'compatible con ancla cronológica'}

def ejecutar(root):
    root=Path(root);F=root/'Fuentes';O=root/'Resultados';O.mkdir(exist_ok=True)
    manifest=json.loads((F/'Manifest_fuentes.json').read_text(encoding='utf-8'))
    for name,digest in manifest['fuentes_sha256'].items():assert hashlib.sha256((F/name).read_bytes()).hexdigest()==digest
    k=pd.read_csv(F/'K_oficiales_tres_PT.csv',dtype={'mis':str},float_precision='round_trip')
    comparisons=pd.read_csv(F/'K_comparacion_tres_PT.csv',dtype={'mis':str},float_precision='round_trip')
    assert k.rol_resultado.eq('oficial').all() and k[k.mis.eq('5c')].T_ka.eq(106).all()
    p=pd.read_csv(F/'perfil_R06_utilizado.csv',float_precision='round_trip')
    assert np.all(np.diff(p.chi_m)>0)
    g=gpd.read_file(F/'cortereloca.shp');assert len(g)==1 and g.geometry.iloc[0].geom_type=='Point' and g.crs is not None
    g=g.to_crs(32718);pt=g.geometry.iloc[0];cut=project(p,[pt.x,pt.y])
    manual=pd.read_csv(F/'kp_manuales_28_09.csv',sep=';',decimal=',');manual=manual[pd.to_numeric(manual.id,errors='coerce').eq(6)]
    candidates=[]
    for v in manual.itertuples():
        q=project(p,[v.x_utm,v.y_utm]);q.update(id_kp=v.id_kp,theta_seleccion_original=v.theta_seleccion);candidates.append(q)
    assert candidates;kp=max(candidates,key=lambda v:v['s_m'])
    rel=k[k.rid.eq(6)];don=k[k.rid.eq(9)]
    assert rel.id_B_corto.eq('N04').all() and len(rel)==3 and len(don)==6
    assert set(don.id_B_corto)=={'A07','A08'} and set(don.mis)=={'5e','5c'}
    rows=[escenario(r,d,cut,kp) for r in rel.itertuples() for d in don.itertuples()];data=pd.DataFrame(rows)
    assert len(data)==18
    curve_frames=[];errors=[]
    for r in rel.itertuples():
        chis=np.unique(np.r_[r.chi_B_m,p.loc[p.chi_m.between(r.chi_B_m,kp['chi_m']),'chi_m'],cut['chi_m'],kp['chi_m']])
        z=np.interp(chis,p.chi_m,p.elevacion_m)
        for d in don.itertuples():
            tc=(cut['chi_m']-r.chi_B_m)/(r.K_m02_anio*FACTOR*1000)
            tau=np.minimum(chis-r.chi_B_m,cut['chi_m']-r.chi_B_m)/(r.K_m02_anio*FACTOR*1000)+np.maximum(chis-cut['chi_m'],0)/(d.K_m02_anio*FACTOR*1000)
            interval_k=np.where((chis[:-1]+chis[1:])/2<=cut['chi_m'],r.K_m02_anio,d.K_m02_anio)
            integral=np.r_[0,np.cumsum(np.diff(chis)/(interval_k*FACTOR*1000))]
            errors.append(float(np.max(np.abs(tau-integral))));assert errors[-1]<1e-8 and tau[0]==0 and (np.diff(tau)>0).all()
            assert abs(tau[chis==cut['chi_m']][0]-tc)<1e-10
            expected=data[(data.id_Reloca==r.id_medicion)&(data.id_donante==d.id_medicion)].iloc[0]
            assert abs(tau[-1]-expected.tau_KP_ka)<1e-8
            curve_frames.append(pd.DataFrame({'id_Reloca':r.id_medicion,'id_donante':d.id_medicion,'MIS_donante':d.mis,
                'chi_m':chis,'elevacion_m':z,'tau_ka':tau,'tau_un_K_Reloca_ka':(chis-r.chi_B_m)/(r.K_m02_anio*FACTOR*1000)}))
    curves=pd.concat(curve_frames,ignore_index=True)
    data.to_csv(O/'Escenarios_dos_K.csv',index=False,encoding='utf-8-sig');data.to_csv(O/'Control_compatibilidad_calibracion.csv',index=False,encoding='utf-8-sig')
    curves.to_csv(O/'Perfiles_chi_altura_tau.csv.gz',index=False,compression='gzip')
    age_rows=[]
    for r in rel.itertuples():
        for d in don.itertuples():
            grid=[(a,a) for a in [-5,0,5]] if d.mis==r.mis else itertools.product([-5,0,5],repeat=2)
            for dr,dd in grid:
                rr=r._replace(T_ka=r.T_ka+dr,K_m02_anio=r.K_m02_anio*r.T_ka/(r.T_ka+dr))
                di=d._replace(T_ka=d.T_ka+dd,K_m02_anio=d.K_m02_anio*d.T_ka/(d.T_ka+dd))
                row=escenario(rr,di,cut,kp);row.update(delta_T_Reloca_ka=dr,delta_T_donante_ka=dd,edad_compartida_MIS=bool(d.mis==r.mis));age_rows.append(row)
    ages=pd.DataFrame(age_rows);assert len(ages)==108
    ages.to_csv(O/'Sensibilidad_edad_y_PT.csv',index=False,encoding='utf-8-sig')
    summary=[]
    for mis,g in data.groupby('MIS_donante'):
        ref=g[g.PT_Reloca.eq('intermedio_real')&g.PT_donante.eq('intermedio_real')].iloc[0]
        summary.append({'MIS_donante':mis,'B_donante':ref.B_donante,'tau_KP_intermedios_ka':ref.tau_KP_ka,
            'tau_corte_intermedio_ka':ref.tau_corte_ka,'tau_corte_KP_intermedio_ka':ref.tau_corte_KP_ka,
            'tau_KP_min_PT_ka':g.tau_KP_ka.min(),'tau_KP_max_PT_ka':g.tau_KP_ka.max(),
            'tau_PT_min_ka':g.tau_PT_dos_K_ka.min(),'tau_PT_max_ka':g.tau_PT_dos_K_ka.max(),
            'n_etiquetas_combinadas':len(g),'n_pares_K_distintos':len(g[['K_Reloca','K_donante']].drop_duplicates()),
            'n_compatibles':int((g.residuo_tau_PT_menos_T_ka.abs()<1e-8).sum()),'n_K_inicial_positivo_posible':int(g.K_inicial_positivo_posible.sum())})
    summary=pd.DataFrame(summary);summary.to_csv(O/'Resumen_tiempos_Reloca.csv',index=False,encoding='utf-8-sig')
    age_summary=ages.groupby('MIS_donante').agg(tau_KP_min_ka=('tau_KP_ka','min'),tau_KP_max_ka=('tau_KP_ka','max'),
            minimo_tiempo_solo_corte_PT_ka=('tau_corte_PT_ka','min'),n_K_inicial_positivo_posible=('K_inicial_positivo_posible','sum')).reset_index()
    age_summary.to_csv(O/'Resumen_sensibilidad_edad.csv',index=False,encoding='utf-8-sig')
    other=comparisons[comparisons.rid.eq(6)]
    extra=pd.DataFrame([escenario(r,d,cut,kp) for r in other.itertuples() for d in don.itertuples()])
    extra.to_csv(O/'Comparacion_otros_B.csv',index=False,encoding='utf-8-sig')
    b=rel.iloc[0];points=[{'punto':'B oficial N04','chi_m':b.chi_B_m,'s_m':b.s_B_m,'z_m':b.z_B_m},
            {'punto':'corte cambio de K',**cut},{'punto':'KP manual final',**kp}]
    for v in rel.itertuples():points.append({'punto':'PT '+v.escenario_PT,'chi_m':v.chi_PT_m,'s_m':v.s_PT_m,'z_m':v.H_m})
    pd.DataFrame(points).to_csv(O/'Puntos_control.csv',index=False,encoding='utf-8-sig')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    colors={'5e':'#9a168d','5c':'#17844e'}
    def plot_chi(ax):
        xx=np.unique(np.r_[b.chi_B_m,p.loc[p.chi_m.between(b.chi_B_m,kp['chi_m']),'chi_m'],cut['chi_m'],kp['chi_m']]);zz=np.interp(xx,p.chi_m,p.elevacion_m)
        for mask,c,label in [(xx<=cut['chi_m'],'#2376a5','K efectivo de Reloca'),(xx>=cut['chi_m'],'#d57920','K transferido de Curanilahue')]:
            ax.plot(xx[mask]-b.chi_B_m,zz[mask],color=c,lw=1.7,label=label)
        for label,x,y,marker,c in [('N04',b.chi_B_m,b.z_B_m,'s','#2376a5'),('Corte',cut['chi_m'],cut['z_m'],'D','#d57920'),('KP',kp['chi_m'],kp['z_m'],'*','#222222')]:
            ax.scatter(x-b.chi_B_m,y,s=75,marker=marker,c=c,zorder=5);ax.annotate(label,(x-b.chi_B_m,y),xytext=(5,9),textcoords='offset points')
        ax.scatter(rel.chi_PT_m-b.chi_B_m,rel.H_m,facecolors='none',edgecolors='#a22342',s=60,zorder=4,label='Tres Pₜ de calibración')
        ax.axhline(b.H_m,color='#a22342',lw=.8,ls=':');ax.set_xlabel('χ desde B=N04 (m)');ax.set_ylabel('Altura DEM (m)');ax.grid(alpha=.2);ax.legend(fontsize=9,loc='upper left')
    def plot_tau(ax):
        cb=b.chi_B_m;chis=np.unique(np.r_[cb,p.loc[p.chi_m.between(cb,kp['chi_m']),'chi_m'],cut['chi_m'],kp['chi_m']]);zz=np.interp(chis,p.chi_m,p.elevacion_m)
        for mis,q in data.groupby('MIS_donante'):
            matrix=np.array([np.minimum(chis-cb,cut['chi_m']-cb)/(v.K_Reloca*FACTOR*1000)+np.maximum(chis-cut['chi_m'],0)/(v.K_donante*FACTOR*1000) for v in q.itertuples()])
            ax.fill_betweenx(zz,matrix.min(axis=0),matrix.max(axis=0),color=colors[mis],alpha=.14,label=f'Donante MIS {mis}: extremos Pₜ')
            for j,v in enumerate(q[q.PT_Reloca.eq('intermedio_real')].drop_duplicates('K_donante').itertuples()):
                tt=np.minimum(chis-cb,cut['chi_m']-cb)/(v.K_Reloca*FACTOR*1000)+np.maximum(chis-cut['chi_m'],0)/(v.K_donante*FACTOR*1000)
                ax.plot(tt,zz,color=colors[mis],lw=1.4,ls=['-','--',':'][j],label=f'{mis}: Pₜ donante {v.PT_donante}')
                ax.scatter(tt[-1],zz[-1],color=colors[mis],s=24)
                ax.scatter(v.tau_PT_dos_K_ka,b.H_m,facecolors='none',edgecolors=colors[mis],s=45,zorder=5)
        r=rel[rel.escenario_PT.eq('intermedio_real')].iloc[0];one=(chis-cb)/(r.K_m02_anio*FACTOR*1000)
        ax.plot(one,zz,color='#555',lw=1,ls='-.',label='Referencia: K uniforme de Reloca')
        tc=(cut['chi_m']-cb)/(r.K_m02_anio*FACTOR*1000);ax.scatter(tc,cut['z_m'],color='#d57920',s=60,marker='D',zorder=5)
        ax.axhline(b.H_m,color='#9d3b47',lw=.7,ls=':')
        ax.scatter(r.T_ka,b.H_m,color='#222222',marker='X',s=65,zorder=6,label='Ancla requerida: 125 ka en Pₜ')
        ax.set_xlabel('τ desde B=N04 (ka)');ax.set_ylabel('Altura DEM (m)');ax.grid(alpha=.2);ax.legend(fontsize=8,loc='lower right')
        ax.text(.02,.98,'Escenarios de transferencia\nNo conservan la calibración τ(Pₜ)=125 ka',transform=ax.transAxes,va='top',fontsize=9,bbox={'fc':'white','ec':'#b77979','alpha':.9})
    for name,fn in [('Reloca_chi_altura_20261020',plot_chi),('Reloca_tau_altura_20261020',plot_tau)]:
        fig,ax=plt.subplots(figsize=(11.5,6.5),layout='constrained');fn(ax);fig.suptitle('Reloca · N04 oficial · dos K por tramos')
        for ext in ['png','svg']:fig.savefig(O/(name+'.'+ext),dpi=180,facecolor='white')
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(17,6.8),sharey=True,layout='constrained');plot_chi(axes[0]);plot_tau(axes[1]);fig.suptitle('Reloca · χ–altura y τ–altura · transferencia condicional de K')
    fig.savefig(O/'Reloca_chi_tau_dos_K_20261020.png',dpi=170,facecolor='white');plt.close(fig)
    text='# Interpretación y alcance\n\n'
    text+='Se aplicó el cambio de K de forma continua en el punto señalado. Los K de Reloca y Curanilahue proceden de las calibraciones oficiales actualizadas; se mantienen separados por MIS y Pₜ. No se promedian coeficientes.\n\n'
    for v in summary.itertuples():text+=f'- Donante {v.MIS_donante} ({v.B_donante}): τ(KP) = {v.tau_KP_intermedios_ka:.2f} ka con ambos Pₜ intermedios; extremos entre Pₜ = {v.tau_KP_min_PT_ka:.2f}–{v.tau_KP_max_PT_ka:.2f} ka. τ(Pₜ) pasa a {v.tau_PT_min_ka:.2f}–{v.tau_PT_max_ka:.2f} ka.\n'
    feasible=int(data.K_inicial_positivo_posible.sum());feasible_age=int(ages.K_inicial_positivo_posible.sum())
    text+=f'\nLos tres Pₜ de Reloca están aguas arriba del corte. En {len(data)-feasible} de {len(data)} combinaciones nominales, el tiempo corte–Pₜ por sí solo supera la edad total de calibración: no queda tiempo positivo para B–corte. Por tanto, ni siquiera reajustando libremente un K inicial positivo se puede mantener τ(Pₜ)=125 ka con esos donantes. El control con rangos de edad conserva {feasible_age} combinaciones factibles de {len(ages)}.\n\n'
    text+='Los tiempos hasta el KP calculados aquí son resultados del escenario solicitado, pero no constituyen una solución compatible con el ancla de terraza bajo esta transferencia. La discrepancia apunta a revisar la representatividad espacial del K donante, la interpretación de Pₜ como ancla temporal o el cambio espacial de K. Este cálculo no permite decidir cuál de esas hipótesis falla.\n\n'
    text+='La posición de cortereloca.shp se utiliza como cambio propuesto por el investigador; aquí no se verifica independientemente que coincida con un contacto litológico. El K de Curanilahue es efectivo para su tramo de calibración: su transferencia como K metamórfico es una hipótesis, no una propiedad medida de la litología. Se conservan las elevaciones originales del DEM, sin suavizar sus oscilaciones. El KP se transfiere por coordenadas desde la selección manual, recalculando χ a θ=0,40.\n'
    (O/'Interpretacion_y_alcance.md').write_text(text,encoding='utf-8')
    controls={'fecha_utc':datetime.now(timezone.utc).isoformat(),'fuentes':manifest,'parametros':{'theta':.4,'A0_m2':1e6,'n':1},
        'n_escenarios_N04':len(data),'n_pares_K_distintos':len(data[['K_Reloca','K_donante']].drop_duplicates()),
        'n_comparacion':len(extra),'n_sensibilidad_edad':len(ages),'edades_MIS5e_compartidas':True,
        'max_error_integracion_independiente_ka':max(errors),'tau_continuo_en_corte':True,'tau_monotono':True,
        'n_K_inicial_positivo_posible':feasible,'n_K_inicial_positivo_posible_con_edad':feasible_age,
        'corte':cut,'KP':kp,'B_oficial':'N04','promedios_K':False,'clasificacion_resultado':'Escenario de transferencia incompatible con ancla cronológica original'}
    (O/'Trazabilidad_y_controles.json').write_text(json.dumps(controls,ensure_ascii=False,indent=2),encoding='utf-8')
    html='<!doctype html><meta charset="utf-8"><title>Reloca dos K</title><style>body{font:16px system-ui;max-width:1400px;margin:35px auto;padding:20px;color:#183247}img{max-width:100%}.nota{padding:16px;background:#fff0ed;border-left:5px solid #ab5148}</style><h1>Reloca · N04 · dos K por tramos</h1><p class="nota">El escenario de transferencia es continuo, pero no conserva la edad de calibración de la terraza. Los tiempos se presentan como escenarios condicionales.</p>'
    short=summary[['MIS_donante','B_donante','tau_KP_intermedios_ka','tau_KP_min_PT_ka','tau_KP_max_PT_ka','tau_PT_min_ka','tau_PT_max_ka']].rename(columns={
        'MIS_donante':'MIS donante','B_donante':'B donante','tau_KP_intermedios_ka':'τ KP intermedios (ka)',
        'tau_KP_min_PT_ka':'τ KP mínimo (ka)','tau_KP_max_PT_ka':'τ KP máximo (ka)','tau_PT_min_ka':'τ Pₜ mínimo (ka)','tau_PT_max_ka':'τ Pₜ máximo (ka)'})
    html+='<style>.tabla{overflow-x:auto}table{font-size:14px;border-collapse:collapse}td,th{padding:8px}</style><div class="tabla">'+short.to_html(index=False,float_format=lambda x:f'{x:.2f}')+'</div><img src="Reloca_chi_altura_20261020.png"><img src="Reloca_tau_altura_20261020.png"><p><a href="Interpretacion_y_alcance.md">Interpretación</a> · <a href="Escenarios_dos_K.csv">Todos los escenarios</a> · <a href="Control_compatibilidad_calibracion.csv">Control de calibración</a></p>'
    (O/'ABRIR_RESULTADOS_RELOCA.html').write_text(html,encoding='utf-8')
    return summary,controls
