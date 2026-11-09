
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
from base_dos_K import project,escenario,FACTOR
ESTILO_FIGURAS={}

def calcular(root):
    root=Path(root);F=root/'Fuentes';O=root/'Resultados';O.mkdir(exist_ok=True)
    manifest=json.loads((F/'Manifest_fuentes.json').read_text(encoding='utf-8'))
    for name,digest in manifest['fuentes_sha256'].items():assert hashlib.sha256((F/name).read_bytes()).hexdigest()==digest
    k=pd.read_csv(F/'K_oficiales_tres_PT.csv',dtype={'mis':str},float_precision='round_trip')
    comp=pd.read_csv(F/'K_comparacion_tres_PT.csv',dtype={'mis':str},float_precision='round_trip')
    p=pd.read_csv(F/'perfil_R06_utilizado.csv',float_precision='round_trip')
    trace=json.loads((F/'Trazabilidad_K.json').read_text(encoding='utf-8'))
    assert trace['parametros']['theta']==.4 and trace['parametros']['A0']==1e6 and trace['parametros']['n']==1
    rel=k[k.rid.eq(6)].copy();don=k[k.rid.eq(9)].copy();other=comp[comp.rid.eq(6)].copy()
    assert len(rel)==3 and rel.id_B_corto.eq('N04').all() and len(don)==6
    assert don[don.mis.eq('5c')].T_ka.eq(106).all() and rel.T_ka.eq(125).all()
    g=gpd.read_file(F/'cortereloca.shp');assert g.crs is not None
    fields=[c for c in g.columns if c.lower()=='id'];assert len(fields)==1
    ids=pd.to_numeric(g[fields[0]],errors='raise');g=g.loc[ids.isin([1,2,3,4])].copy();ids=ids.loc[g.index];assert set(ids)=={1,2,3,4} and ids.is_unique and len(g)==4
    g=g.to_crs(32718);cuts={}
    for v,geom in zip(ids,g.geometry):
        assert geom.geom_type=='Point';cuts[int(v)]=project(p,[geom.x,geom.y])
    cuts=dict(sorted(cuts.items()))
    manual=pd.read_csv(F/'kp_manuales_28_09.csv',sep=';',decimal=',');manual=manual[pd.to_numeric(manual.id,errors='coerce').eq(6)]
    kp_list=[]
    for v in manual.itertuples():
        a=project(p,[v.x_utm,v.y_utm]);a.update(id_kp=v.id_kp,theta_seleccion_original=float(v.theta_seleccion));kp_list.append(a)
    assert kp_list;kp=max(kp_list,key=lambda a:a['s_m'])
    assert np.all(np.diff(p.chi_m)>0)
    for cut in cuts.values():assert float(rel.chi_B_m.max())<cut['chi_m']<kp['chi_m']
    rows=[];age_rows=[];comparison=[];curve_frames=[];errors=[]
    for cid,cut in cuts.items():
        for r in rel.itertuples():
            chis=np.unique(np.r_[r.chi_B_m,p.loc[p.chi_m.between(r.chi_B_m,kp['chi_m']),'chi_m'],cut['chi_m'],kp['chi_m']])
            z=np.interp(chis,p.chi_m,p.elevacion_m)
            for d in don.itertuples():
                row=escenario(r,d,cut,kp)|{'id_corte':cid};rows.append(row)
                tau=np.minimum(chis-r.chi_B_m,cut['chi_m']-r.chi_B_m)/(r.K_m02_anio*FACTOR*1000)+np.maximum(chis-cut['chi_m'],0)/(d.K_m02_anio*FACTOR*1000)
                interval_k=np.where((chis[:-1]+chis[1:])/2<=cut['chi_m'],r.K_m02_anio,d.K_m02_anio)
                independent=np.r_[0,np.cumsum(np.diff(chis)/(interval_k*FACTOR*1000))]
                error=float(np.max(np.abs(tau-independent)));errors.append(error)
                assert error<1e-8 and tau[0]==0 and (np.diff(tau)>0).all()
                assert abs(tau[-1]-row['tau_KP_ka'])<1e-8
                assert abs(tau[chis==cut['chi_m']][0]-row['tau_corte_ka'])<1e-10
                curve_frames.append(pd.DataFrame({'id_corte':cid,'id_Reloca':r.id_medicion,'id_donante':d.id_medicion,
                    'MIS_donante':d.mis,'chi_m':chis,'elevacion_m':z,'tau_ka':tau,
                    'tau_un_K_Reloca_ka':(chis-r.chi_B_m)/(r.K_m02_anio*FACTOR*1000)}))
                grid=[(a,a) for a in [-5,0,5]] if r.mis==d.mis else itertools.product([-5,0,5],repeat=2)
                for dr,dd in grid:
                    rr=r._replace(T_ka=r.T_ka+dr,K_m02_anio=r.K_m02_anio*r.T_ka/(r.T_ka+dr))
                    di=d._replace(T_ka=d.T_ka+dd,K_m02_anio=d.K_m02_anio*d.T_ka/(d.T_ka+dd))
                    age_rows.append(escenario(rr,di,cut,kp)|{'id_corte':cid,'delta_T_Reloca_ka':dr,'delta_T_donante_ka':dd,'edad_compartida_MIS':r.mis==d.mis})
        for r in other.itertuples():
            for d in don.itertuples():comparison.append(escenario(r,d,cut,kp)|{'id_corte':cid})
    data=pd.DataFrame(rows);ages=pd.DataFrame(age_rows);curves=pd.concat(curve_frames,ignore_index=True)
    assert len(data)==72 and len(ages)==432
    assert data.groupby('id_corte').size().eq(18).all()
    assert not data.duplicated(['id_corte','id_Reloca','id_donante']).any()
    data['conserva_calibracion_original']=data.residuo_tau_PT_menos_T_ka.abs()<1e-8
    for cid,cut in cuts.items():
        lower=data[data.id_corte.eq(cid)&~data.PT_aguas_arriba_del_corte]
        assert lower.conserva_calibracion_original.all()
    data.to_csv(O/'Escenarios_72_dos_K_cuatro_cortes.csv',index=False,encoding='utf-8-sig')
    ages.to_csv(O/'Sensibilidad_432_edad_y_PT.csv',index=False,encoding='utf-8-sig')
    curves.to_csv(O/'Perfiles_chi_altura_tau_cuatro_cortes.csv.gz',index=False,compression='gzip')
    pd.DataFrame(comparison).to_csv(O/'Comparacion_B_no_oficiales.csv',index=False,encoding='utf-8-sig')
    summary=[]
    for (cid,mis),q in data.groupby(['id_corte','MIS_donante']):
        ref=q[q.PT_Reloca.eq('intermedio_real')&q.PT_donante.eq('intermedio_real')].iloc[0]
        cut=cuts[int(cid)];r=rel[rel.escenario_PT.eq('intermedio_real')].iloc[0]
        status='los tres Pₜ conservan T' if q.conserva_calibracion_original.all() else 'algunos Pₜ conservan T' if q.conserva_calibracion_original.any() else 'ningún Pₜ conserva T'
        summary.append({'id_corte':int(cid),'MIS_donante':mis,'B_donante':ref.B_donante,'altura_corte_m':cut['z_m'],
            'distancia_B_corte_km':(cut['s_m']-r.s_B_m)/1000,'chi_corte_m':cut['chi_m'],
            'tau_hasta_corte_intermedio_ka':ref.tau_corte_ka,
            'tau_hasta_corte_min_PT_ka':q.tau_corte_ka.min(),'tau_hasta_corte_max_PT_ka':q.tau_corte_ka.max(),
            'tau_corte_KP_intermedio_ka':ref.tau_corte_KP_ka,'tau_KP_intermedios_ka':ref.tau_KP_ka,
            'tau_KP_min_PT_ka':q.tau_KP_ka.min(),'tau_KP_max_PT_ka':q.tau_KP_ka.max(),
            'tau_PT_min_ka':q.tau_PT_dos_K_ka.min(),'tau_PT_max_ka':q.tau_PT_dos_K_ka.max(),
            'residuo_PT_min_ka':q.residuo_tau_PT_menos_T_ka.min(),'residuo_PT_max_ka':q.residuo_tau_PT_menos_T_ka.max(),
            'n_PT_Reloca_preservados':q[q.conserva_calibracion_original].PT_Reloca.nunique(),
            'n_combinaciones_compatibles':int(q.conserva_calibracion_original.sum()),'n_combinaciones':len(q),
            'n_K_inicial_positivo_posible':int(q.K_inicial_positivo_posible.sum()),'estado_calibracion':status})
    summary=pd.DataFrame(summary);summary.to_csv(O/'Resumen_cuatro_cortes.csv',index=False,encoding='utf-8-sig')
    projection=pd.DataFrame([{'id_corte':cid,**cut} for cid,cut in cuts.items()]);projection.to_csv(O/'Cortes_proyectados.csv',index=False,encoding='utf-8-sig')
    age_summary=ages.groupby(['id_corte','MIS_donante']).agg(tau_KP_min_ka=('tau_KP_ka','min'),tau_KP_max_ka=('tau_KP_ka','max'),
        n_K_inicial_positivo_posible=('K_inicial_positivo_posible','sum')).reset_index()
    age_summary.to_csv(O/'Resumen_sensibilidad_edad.csv',index=False,encoding='utf-8-sig')
    plotear(root,p,rel,don,cuts,kp,data,summary,estilo=ESTILO_FIGURAS)
    report={'fecha_utc':datetime.now(timezone.utc).isoformat(),'fuentes':manifest,'parametros':{'theta':.4,'A0_m2':1e6,'n':1},
        'B_oficial':'N04','n_cortes':4,'ids_cortes':list(cuts),'cortes':cuts,'KP':kp,
        'n_escenarios_nominales':len(data),'n_sensibilidad_edad':len(ages),'n_comparacion':len(comparison),
        'max_error_integracion_independiente_ka':max(errors),'tau_continuo_en_cada_corte':True,'promedios_K':False,
        'cronologia_publicada':{'5e':[125,5],'5c':[106,5]},'edades_dentro_de_MIS5e_compartidas':True,
        'compatibles_por_corte_MIS':summary[['id_corte','MIS_donante','n_combinaciones_compatibles','n_PT_Reloca_preservados']].to_dict('records')}
    (O/'Trazabilidad_y_controles.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    text='# Reloca: cuatro cambios propuestos de K\n\nB oficial N04. Hasta cada corte se usa K calibrado de Reloca; desde ese corte hasta el KP manual se usa cada K de Curanilahue por separado. Los cortes son cuatro alternativas independientes, no cuatro cambios consecutivos dentro del mismo modelo. Sin promedios de MIS ni de Pₜ.\n\n'
    for cid in cuts:
        v=summary[summary.id_corte.eq(cid)]
        text+=f'## Corte {cid}\n\n'
        for row in v.itertuples():text+=f'- Donante {row.MIS_donante}: τ hasta corte = {row.tau_hasta_corte_intermedio_ka:.2f} ka; corte–KP = {row.tau_corte_KP_intermedio_ka:.2f} ka; total KP = {row.tau_KP_intermedios_ka:.2f} ka. Extremos de Pₜ: {row.tau_KP_min_PT_ka:.2f}–{row.tau_KP_max_PT_ka:.2f} ka. Calibración: {row.estado_calibracion}.\n'
        text+='\n'
    text+='La conservación de τ(Pₜ)=T es un control algebraico, no una validación geológica independiente de la posición del corte o del K transferido. Un corte situado aguas arriba de los tres Pₜ conserva automáticamente el ancla, porque el cambio de K ocurre después del tramo calibrado. Si el corte cae entre los Pₜ, cada escenario se verifica por separado.\n\n'
    text+='Los perfiles conservan las elevaciones del DEM. La selección manual del KP se transfiere por sus coordenadas; χ se recalcula con θ=0,40. El punto original fue seleccionado con θ=0,45. Las cotas de los cortes son las del perfil interpolado y no de una altura digitada. Los rangos cronológicos se comparten entre calibraciones del mismo MIS y se muestran aparte de los rangos por Pₜ. No son intervalos de confianza.\n'
    (O/'Interpretacion_y_alcance.md').write_text(text,encoding='utf-8')
    html=['<!doctype html><html lang="es"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Reloca · cuatro cortes</title><style>body{font:16px system-ui;color:#183247;max-width:1500px;margin:30px auto;padding:15px}img{max-width:100%;height:auto}.par{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:15px}.tabla{overflow-x:auto}td,th{padding:8px;text-align:right}table{font-size:14px;border-collapse:collapse}section{margin:35px 0;padding-top:15px;border-top:1px solid #ddd}nav a{margin-right:18px}a{color:#1764a0}.nota{padding:15px;background:#eef4f8}@media(max-width:850px){.par{grid-template-columns:1fr}}</style><h1>Reloca · cuatro cortes alternativos</h1><p class="nota">B oficial N04. Cada ejemplo usa K de Reloca hasta su corte y K de Curanilahue desde allí al KP. MIS 5e y 5c y los tres Pₜ se mantienen separados. Son cuatro alternativas independientes.</p><nav>'+''.join(f'<a href="#corte{cid}">Corte {cid}</a>' for cid in cuts)+'</nav>']
    short=tabla_resumen(summary)
    html+=['<h2>Comparación de tiempos (ka)</h2><div class="tabla">'+short.to_html(index=False,float_format=lambda x:f'{x:.2f}')+'</div>','<img src="Reloca_cuatro_cortes_chi_altura_20261020.png">']
    for cid in cuts:
        status=summary[summary.id_corte.eq(cid)].estado_calibracion.iloc[0]
        html.append(f'<section id="corte{cid}"><h2>Corte {cid}</h2><p>Control de calibración: {status}. Las líneas usan Pₜ intermedio de Reloca; las líneas distinguen los Pₜ con edades nominales, sin rellenos.</p><div class="par"><a href="Corte_{cid}_chi_altura.png"><img src="Corte_{cid}_chi_altura.png"></a><a href="Corte_{cid}_tau_altura.png"><img src="Corte_{cid}_tau_altura.png"></a></div></section>')
    html+=['<p>La conservación del ancla de terraza no verifica de forma independiente la transferencia del K donante ni un contacto litológico.</p><p><a href="Resumen_cuatro_cortes.csv">Tabla completa</a> · <a href="Escenarios_72_dos_K_cuatro_cortes.csv">72 escenarios nominales</a> · <a href="Interpretacion_y_alcance.md">Interpretación y alcance</a></p></html>']
    (O/'ABRIR_RELOCA_CUATRO_CORTES.html').write_text('\n'.join(html),encoding='utf-8')
    return summary,report

def tabla_resumen(summary):
    return summary[['id_corte','MIS_donante','tau_hasta_corte_intermedio_ka','tau_corte_KP_intermedio_ka','tau_KP_intermedios_ka','tau_KP_min_PT_ka','tau_KP_max_PT_ka','n_PT_Reloca_preservados']].rename(columns={
        'id_corte':'Corte','MIS_donante':'MIS donante','tau_hasta_corte_intermedio_ka':'τ hasta corte (ka)',
        'tau_corte_KP_intermedio_ka':'τ corte–KP (ka)','tau_KP_intermedios_ka':'τ total KP (ka)',
        'tau_KP_min_PT_ka':'τ total mínimo Pₜ','tau_KP_max_PT_ka':'τ total máximo Pₜ','n_PT_Reloca_preservados':'Pₜ conservados de 3'})

def plotear(root,p,rel,don,cuts,kp,data,summary,estilo=None):
    O=root/'Resultados';estilo=estilo or {};colors=estilo.get('COLORES_MIS',{'5e':'#9a168d','5c':'#17844e'});cutcolors=estilo.get('COLORES_CORTES',{1:'#d18315',2:'#477bc1',3:'#b13d53',4:'#6c4aa2'})
    ref=rel[rel.escenario_PT.eq('intermedio_real')].iloc[0];cb=ref.chi_B_m;ck=kp['chi_m']
    plt.rcParams.update({'font.family':estilo.get('FUENTE','DejaVu Sans'),'font.size':estilo.get('TAMANO_FUENTE',10),'axes.spines.top':False,'axes.spines.right':False})
    def chi_plot(ax,cid=None):
        chis=np.unique(np.r_[cb,p.loc[p.chi_m.between(cb,ck),'chi_m'],*[c['chi_m'] for c in cuts.values()],ck]);zz=np.interp(chis,p.chi_m,p.elevacion_m)
        if cid is None:ax.plot(chis-cb,zz,color='#647988',lw=1.4,label='Perfil DEM')
        else:
            ct=cuts[cid]['chi_m']
            for mask,col,label in [(chis<=ct,'#2376a5','K de Reloca'),(chis>=ct,'#d57920','K de Curanilahue')]:ax.plot(chis[mask]-cb,zz[mask],c=col,lw=1.6,label=label)
        for number,cut in cuts.items():
            if cid is not None and cid!=number:continue
            ax.scatter(cut['chi_m']-cb,cut['z_m'],c=cutcolors[number],s=55,marker=estilo.get('SIMBOLO_CORTE','D'),zorder=5)
            ax.annotate(f'Corte {number}',(cut['chi_m']-cb,cut['z_m']),xytext=(7,10),textcoords='offset points',color=cutcolors[number])
        ax.scatter(rel.chi_PT_m-cb,rel.H_m,facecolors='none',edgecolors='#c12b49',s=55,zorder=6,label='Tres Pₜ de calibración')
        ax.axhline(ref.H_m,color='#c12b49',ls=':',lw=.7)
        for x,y,label,marker in [(0,ref.z_B_m,'N04',estilo.get('SIMBOLOS',{}).get('B','s')),(ck-cb,kp['z_m'],'KP',estilo.get('SIMBOLOS',{}).get('KP','*'))]:
            ax.scatter(x,y,c='#222',s=55,marker=marker,zorder=5);ax.annotate(label,(x,y),xytext=(5,7),textcoords='offset points')
        ax.set_xlabel('χ desde B=N04 (m)');ax.set_ylabel('Altura DEM (m)');ax.grid(alpha=.2);ax.legend(loc='upper left',fontsize=8)
        ax.set_ylim(bottom=0)
    def tau_plot(ax,cid,shared_x=None):
        cut=cuts[cid];ct=cut['chi_m'];chis=np.unique(np.r_[cb,p.loc[p.chi_m.between(cb,ck),'chi_m'],ct,ck]);zz=np.interp(chis,p.chi_m,p.elevacion_m)
        for mis,q in data[data.id_corte.eq(cid)].groupby('MIS_donante'):
            matrix=np.array([np.minimum(chis-cb,ct-cb)/(v.K_Reloca*FACTOR*1000)+np.maximum(chis-ct,0)/(v.K_donante*FACTOR*1000) for v in q.itertuples()])
            for j,v in enumerate(q[q.PT_Reloca.eq('intermedio_real')].drop_duplicates('K_donante').itertuples()):
                tt=np.minimum(chis-cb,ct-cb)/(v.K_Reloca*FACTOR*1000)+np.maximum(chis-ct,0)/(v.K_donante*FACTOR*1000)
                label=f'{mis} · donante {v.PT_donante}' if len(q.K_donante.unique())>1 else f'{mis} · tres Pₜ coincidentes'
                ax.plot(tt,zz,color=colors[mis],lw=1.3,ls=['-','--',':'][j],label=label)
                ax.scatter(tt[-1],zz[-1],s=20,color=colors[mis]);ax.scatter(v.tau_PT_dos_K_ka,ref.H_m,facecolors='none',edgecolors=colors[mis],s=40,zorder=5)
        one=(chis-cb)/(ref.K_m02_anio*FACTOR*1000);ax.plot(one,zz,c='#777',ls='-.',lw=.9,label='Referencia: K uniforme Reloca')
        tc=(ct-cb)/(ref.K_m02_anio*FACTOR*1000);ax.scatter(tc,cut['z_m'],c=cutcolors[cid],marker=estilo.get('SIMBOLO_CORTE','D'),s=55,zorder=5)
        ax.scatter(ref.T_ka,ref.H_m,c='#222',marker='X',s=60,zorder=6,label='Ancla de referencia: 125 ka')
        ax.axhline(ref.H_m,c='#c12b49',ls=':',lw=.7)
        ax.set_xlabel('τ desde B=N04 (ka)');ax.set_ylabel('Altura DEM (m)');ax.grid(alpha=.2);ax.set_xlim(left=0);ax.set_ylim(bottom=0)
        if shared_x:ax.set_xlim(0,shared_x)
        ax.legend(loc='lower right',fontsize=7)
        state=summary[summary.id_corte.eq(cid)].estado_calibracion.iloc[0]
        ax.set_title(f'Corte {cid} · {state}',fontsize=11)
    fig,ax=plt.subplots(figsize=(11,6),layout='constrained');chi_plot(ax);fig.suptitle('Reloca · ubicación de los cuatro cortes alternativos')
    fig.savefig(O/'Reloca_cuatro_cortes_chi_altura_20261020.png',dpi=estilo.get('DPI',180),facecolor='white');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(15,11.5),sharex=True,sharey=True,layout='constrained')
    for ax,cid in zip(axes.flat,cuts):tau_plot(ax,cid,shared_x=float(data.tau_KP_ka.max())*1.04)
    fig.suptitle('Reloca · τ–altura para cuatro cambios alternativos de K\nN04 oficial · K donante 5e y 5c separados · edades nominales')
    fig.savefig(O/'Reloca_cuatro_cortes_tau_altura_20261020.png',dpi=estilo.get('DPI',160),facecolor='white');plt.close(fig)
    for cid in cuts:
        for name,fn in [('chi_altura',chi_plot),('tau_altura',tau_plot)]:
            fig,ax=plt.subplots(figsize=(10,6),layout='constrained');fn(ax,cid)
            fig.suptitle(f'Reloca · corte {cid} · N04 oficial')
            for ext in estilo.get('FORMATOS',['png','svg']):fig.savefig(O/f'Corte_{cid}_{name}.{ext}',dpi=estilo.get('DPI',165),facecolor='white')
            plt.close(fig)
