"""Reproduce the thesis synthesis from the explicit 22 October derived-data snapshot.

This is not a rerun of DEM routing or an independent validation of tectonic ages.
K_G denotes propagation calibration; K_P denotes slope/rate calibration. Units:
K: m**0.2/year; chi: m with A0=1e6 m2; tau: ka; modeled incision: mm/year.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
THETA=0.40
A0=1_000_000.0
FACTOR=A0**THETA

def verificar_archivos():
    manifest=json.loads((ROOT/'datos_sintesis/manifest.json').read_text(encoding='utf8'))
    for name,digest in manifest['sha256'].items():
        path=ROOT/'datos_sintesis'/name
        if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            raise ValueError(f'Cambió el insumo {name}; revisar procedencia antes de calcular.')
    return manifest

def cargar():
    verificar_archivos()
    read=lambda name:pd.read_csv(ROOT/'datos_sintesis'/name,float_precision='round_trip')
    cal=read('Calibraciones.csv')
    assert len(cal)==24 and cal.groupby(['rid','mis']).chi_B_m.nunique().eq(1).all()
    assert cal.groupby('mis').T_ka.first().to_dict()=={'5c':106.0,'5e':125.0}
    # Three labels can refer to one physical crossing, so each physical point counts once.
    cal['ptkey']=cal.chi_PT_m.round(7)
    physical=cal.drop_duplicates(['rid','mis','ptkey']).copy()
    assert len(physical)==18
    means=physical.groupby(['rid','mis'],as_index=False).agg(
        rio=('rio','first'),K_G=('K_George','mean'),K_P=('K_pendiente','mean'),
        latitud_sur=('latitud_sur','first'),n_PT_fisicos=('ptkey','size'))
    kp=read('KP.csv').set_index('rid')
    origins=read('Geometria_origenes.csv')
    profiles={rid:read(f'R{rid:02d}.csv') for rid in range(1,19)}
    for rid,p in profiles.items():
        assert np.isfinite(p.to_numpy(float)).all(),rid
        assert (np.diff(p.distancia_m)>0).all() and p.area_m2.gt(0).all(),rid
        weights=p.area_m2.to_numpy()**-THETA
        chi=FACTOR*np.r_[0.0,np.cumsum(0.5*(weights[:-1]+weights[1:])*np.diff(p.distancia_m))]
        assert np.allclose(chi,p.chi_m-p.chi_m.iloc[0],rtol=1e-10,atol=1e-6),rid
    return cal,physical,means,kp,origins,profiles,read('Cortes_Reloca.csv').set_index('id_corte')

def calcular():
    cal,physical,means,kp,origins,profiles,cuts=cargar()
    rows=[]
    for o in origins.itertuples():
        rid=int(o.rid)
        if rid not in kp.index or rid==1:
            continue  # R01 has no assigned K; R07 has no accepted manual KP.
        donor=4 if rid in (2,3) else 9 if rid>=10 else rid
        # Reloca always keeps its locally calibrated 5e lower domain.
        mis='5e' if rid==6 else o.mis
        q=means[(means.rid==donor)&(means.mis==mis)]
        assert len(q)==1,(rid,mis)
        row=q.iloc[0]
        result={'rid':rid,'rio':o.rio,'mis':o.mis,'donante_K':donor,
                'lat_KP_sur':o.lat_KP_sur,'z_KP_m':kp.loc[rid,'z_m'],
                'chi_origen_m':o.chi_origen_m,'K_G':row.K_G,'K_P':row.K_P}
        for method in ('G','P'):
            k=float(row['K_'+method])
            if rid==6:
                cut=cuts.loc[2]
                upper=means[(means.rid==9)&(means.mis==o.mis)].iloc[0]['K_'+method]
                result['K_superior_'+method]=upper
                result['tau_'+method+'_ka']=(cut.chi_m-o.chi_origen_m)/(k*FACTOR*1000)+(kp.loc[rid,'chi_m']-cut.chi_m)/(upper*FACTOR*1000)
            else:
                result['tau_'+method+'_ka']=(kp.loc[rid,'chi_m']-o.chi_origen_m)/(k*FACTOR*1000)
        rows.append(result)
    times=pd.DataFrame(rows)
    south=[]
    for mis in ('5c','5e'):
        q=times[(times.rid>=10)&(times.mis==mis)]
        assert len(q)==9
        for method in ('G','P'):
            values=q['tau_'+method+'_ka']
            south.append({'mis':mis,'metodo':method,'n':len(values),'min_ka':values.min(),
                          'max_ka':values.max(),'media_ka':values.mean(),'mediana_ka':values.median()})
    # The C2 model may not preserve all original calibration anchors. Report its residuals.
    residuals=[]
    lower=physical[physical.rid==6]
    cut=cuts.loc[2]
    for mis in ('5c','5e'):
        upper=means[(means.rid==9)&(means.mis==mis)].K_G.iloc[0]
        for row in lower.itertuples():
            age=(min(row.chi_PT_m,cut.chi_m)-row.chi_B_m)/(row.K_George*FACTOR*1000)+max(row.chi_PT_m-cut.chi_m,0)/(upper*FACTOR*1000)
            residuals.append({'mis_superior':mis,'s_PT_m':row.s_PT_m,'tau_PT_ka':age,
                              'referencia_ka':row.T_ka,'residuo_ka':age-row.T_ka})
    rates=[]
    for row in times.itertuples():
        p=profiles[row.rid];point=kp.loc[row.rid]
        q=p[p.distancia_m.between(point.s_m-2000,point.s_m)]
        if len(q)<5:
            continue
        fit=np.polyfit(q.chi_m,q.elevacion_m,1)
        predicted=np.polyval(fit,q.chi_m)
        r2=1-((q.elevacion_m-predicted)**2).sum()/((q.elevacion_m-q.elevacion_m.mean())**2).sum()
        ksn=fit[0]*FACTOR
        kg=row.K_superior_G if row.rid==6 else row.K_G
        rates.append({'rid':row.rid,'rio':row.rio,'mis':row.mis,'latitud_sur':row.lat_KP_sur,
                      'K_G':kg,'ksn_KP_m08':ksn,'R2':r2,'n':len(q),
                      'ventana_m':q.distancia_m.max()-q.distancia_m.min(),
                      'E_mm_anio':1000*kg*ksn})
    rates=pd.DataFrame(rates)
    assert rates.E_mm_anio.gt(0).all()
    return {'K_promedio':means,'Tiempos_C2':times,'Resumen_sur':pd.DataFrame(south),
            'Residuos_PT_C2':pd.DataFrame(residuals),'Tasas_bajo_KP':rates},profiles,kp

def verificar_referencia(tables):
    ref=pd.read_csv(ROOT/'datos_sintesis/Referencia_tiempos.csv',float_precision='round_trip')
    got=tables['Tiempos_C2'].set_index(['rid','mis'])
    for r in ref.itertuples():
        if not np.isfinite(r.tau_George_ka):continue
        row=got.loc[(r.rid,r.mis)]
        assert np.isclose(row.tau_G_ka,r.tau_George_ka,rtol=1e-10,atol=1e-7),(r.rid,r.mis,'G')
        assert np.isclose(row.tau_P_ka,r.tau_pendiente_ka,rtol=1e-10,atol=1e-7),(r.rid,r.mis,'P')
    reference=pd.read_csv(ROOT/'datos_sintesis/Referencia_tasas.csv',float_precision='round_trip').set_index(['rid','mis'])
    got=tables['Tasas_bajo_KP'].set_index(['rid','mis'])
    assert set(got.index)==set(reference.index)
    assert np.allclose(got.E_mm_anio,reference.loc[got.index].E_mm_anio,rtol=1e-10,atol=1e-10)

def graficar(tables,profiles,kp,output):
    colors={'5c':'#2471a3','5e':'#c05b26'}
    paths=[]
    def save(fig,name):
        for ext in ('png','svg'):
            p=output/(name+'.'+ext);fig.savefig(p,dpi=200,bbox_inches='tight');paths.append(p)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),sharey=True)
    for ax,method in zip(axes,('G','P')):
        for mis in ('5c','5e'):
            q=tables['K_promedio'].query('mis == @mis')
            ax.scatter(q.latitud_sur,q['K_'+method],color=colors[mis],label=mis)
            if len(q)>2:
                fit=np.polyfit(q.latitud_sur,q['K_'+method],1);x=np.linspace(q.latitud_sur.min(),q.latitud_sur.max())
                ax.plot(x,np.polyval(fit,x),color=colors[mis],ls='--',lw=1)
        ax.set(xlabel='Latitud del río (°S)',ylabel=f'$K_{method}$ (m⁰·²/año)');ax.legend();ax.grid(alpha=.15)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    fig.tight_layout();save(fig,'K_latitud')
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True)
    for i,mis in enumerate(('5c','5e')):
        q=tables['Tiempos_C2'].query('mis == @mis')
        for j,method in enumerate(('G','P')):
            ax=axes[i,j];ax.scatter(q.lat_KP_sur,q['tau_'+method+'_ka'],color=colors[mis])
            for r in q.itertuples():ax.annotate(f'{r.rid:02d}',(r.lat_KP_sur,getattr(r,'tau_'+method+'_ka')),fontsize=7,xytext=(3,4),textcoords='offset points')
            ax.set(title=f'MIS {mis} · $K_{method}$',ylabel='τ (ka)',xlabel='Latitud del KP (°S)');ax.grid(alpha=.15);ax.margins(y=.12)
    fig.tight_layout();save(fig,'Tau_MIS_metodo_C2')
    fig,axes=plt.subplots(1,2,figsize=(10,4.5),sharey=True)
    for ax,mis in zip(axes,('5c','5e')):
        q=tables['Tasas_bajo_KP'].query('mis == @mis')
        ax.scatter(q.latitud_sur,q.E_mm_anio,color=colors[mis])
        for r in q.itertuples():
            offset={17:(-16,-12),18:(4,10),12:(-15,-12)}.get(r.rid,(3,5))
            ax.annotate(f'{r.rid:02d}',(r.latitud_sur,r.E_mm_anio),fontsize=7,xytext=offset,textcoords='offset points')
        ax.set(title=f'$K_G$ MIS {mis}',xlabel='Latitud del KP (°S)');ax.grid(alpha=.15);ax.margins(x=.1,y=.12)
    axes[0].set_ylabel('Incisión modelada $K_G k_{sn}$ (mm/año)\nU aparente sólo bajo equilibrio local')
    fig.tight_layout();save(fig,'Tasas_KG_sin_referencias_costeras')
    fig=plt.figure(figsize=(11,7.5));ax=fig.add_subplot(111,projection='3d')
    palette=plt.get_cmap('tab20')
    for rid,p in profiles.items():
        q=p.iloc[::max(1,len(p)//500)];color=palette((rid-1)%20)
        ax.plot((q.x_utm-700000)/1000,(q.y_utm-5900000)/1000,q.elevacion_m,c=color,label=f'R{rid:02d}')
        if rid in kp.index:
            v=kp.loc[rid];ax.scatter((v.x_perfil_m-700000)/1000,(v.y_perfil_m-5900000)/1000,v.z_m,marker='v',c=[color],s=35)
    ax.view_init(elev=28,azim=260);ax.set_box_aspect((1,2.1,.85))
    ax.set_xlabel('Este UTM −700 km (km)',fontsize=9,labelpad=8)
    ax.set_ylabel('Norte UTM −5.900 km (km)',fontsize=9,labelpad=10)
    ax.set_zlabel('Elevación DEM (m)',fontsize=9,labelpad=7)
    ax.tick_params(labelsize=8)
    ax.legend(loc='center left',bbox_to_anchor=(1.05,.5),fontsize=7)
    fig.text(.05,.02,'Vista hacia N10°E · EPSG:32718 · Escala vertical exagerada · Sin hillshade',fontsize=8)
    save(fig,'Perfiles_3D_N10E_sin_hillshade')
    return paths

def ejecutar(output=None):
    if output is None:
        stamp=datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%fZ')
        output=ROOT/'Productos/sintesis_20261022'/stamp
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    tables,profiles,kp=calcular();verificar_referencia(tables)
    for name,table in tables.items():table.to_csv(output/(name+'.csv'),index=False,encoding='utf8')
    figures=graficar(tables,profiles,kp,output)
    report={'theta':THETA,'A0_m2':A0,'n':1,'reloca':'C2; 5e lower domain; 5c or 5e upper donor',
            'parity_with_saved_thesis_summary':True,'full_dem_rerun':False,
            'external_unpublished_data_included':False,'hillshade':False,
            'figures':[p.name for p in figures],
            'uncertainty':'Scenario ranges are not confidence intervals. U=E requires local equilibrium.'}
    (output/'verificacion.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')
    return output,tables

if __name__=='__main__':
    out,tables=ejecutar()
    print(tables['Resumen_sur'].to_string(index=False))
    print('Resultados nuevos:',out)
