
import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import LineString,Point

GRANT=pd.DataFrame(columns=['muestra','K','error_reportado'])
COLORS=['#1764a0','#e17c13','#7e42a7','#158663','#c64359','#82602b']

def limites_grant():
    if GRANT.empty:
        return dict.fromkeys(['central_min','central_max','outer_min','outer_max'], np.nan)
    return {'central_min':float(GRANT.K.min()),'central_max':float(GRANT.K.max()),
            'outer_min':float((GRANT.K-GRANT.error_reportado).min()),
            'outer_max':float((GRANT.K+GRANT.error_reportado).max())}

def clasificar_grant(datos):
    if GRANT.empty:
        d=datos.copy()
        d['comparacion_Grant_central']='sin referencia externa distribuida'
        d['dentro_Grant_con_errores']=pd.NA
        return d
    d=datos.copy();b=limites_grant()
    d['comparacion_Grant_central']=np.where(d.K_m02_anio<b['central_min'],'debajo',np.where(d.K_m02_anio>b['central_max'],'encima','dentro'))
    d['dentro_Grant_con_errores']=d.K_m02_anio.between(b['outer_min'],b['outer_max'])
    return d

def calcular_tau(resultados,perfiles,kp_path,folder,A0=1e6,theta=.4,mostrar_todas=True):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    kp=pd.read_csv(kp_path,sep=';',decimal=',')
    curvas=[];envolventes=[];valoreskp=[];figures=[];verificaciones=[]
    for rid,group in resultados.groupby('rid',sort=True):
        p=perfiles[int(rid)];line=LineString(p[['x_utm','y_utm']].to_numpy())
        geom_s=np.r_[0,np.cumsum(np.hypot(np.diff(p.x_utm),np.diff(p.y_utm)))]
        manual=kp[pd.to_numeric(kp.id,errors='coerce').eq(int(rid))].copy()
        kpinfo=[]
        for _,k in manual.iterrows():
            point=Point(k.x_utm,k.y_utm);s=line.project(point);distance=point.distance(line)
            chi=float(np.interp(s,geom_s,p.chi_m));z=float(np.interp(s,geom_s,p.elevacion_m))
            kpinfo.append({'id_kp':k.id_kp,'chi_kp_m':chi,'z_kp_m':z,'separacion_cauce_m':distance,'valido':distance<=35})
        fig,(az,at)=plt.subplots(1,2,figsize=(14,6),layout='constrained')
        for color,(case,q) in zip(COLORS,group.groupby('caso_B',sort=False)):
            b=q.iloc[0];tail=p[p.chi_m>b.chi_B_m+1e-8]
            x=np.r_[b.chi_B_m,tail.chi_m.to_numpy()];z=np.r_[b.z_B_m,tail.elevacion_m.to_numpy()]
            matrix=[]
            for _,v in q.iterrows():
                tau=(x-v.chi_B_m)/(v.K_m02_anio*A0**theta)/1000
                assert tau[0]==0 and np.all(tau>=0)
                anchor=(v.chi_PT_m-v.chi_B_m)/(v.K_m02_anio*A0**theta)/1000
                verificaciones.append(abs(anchor-v.T_ka))
                matrix.append(tau)
                curvas.append(pd.DataFrame({'rid':rid,'id_medicion':v.id_medicion,'caso_B':case,
                    'K_m02_anio':v.K_m02_anio,'chi_m':x,'elevacion_m':z,'tau_ka':tau}))
                if mostrar_todas:
                    az.plot(tau,z,color=color,alpha=.20,lw=.65)
                    at.plot(x,tau,color=color,alpha=.12,lw=.5)
                for k in kpinfo:
                    usable=k['valido'] and k['chi_kp_m']>=v.chi_B_m
                    t=(k['chi_kp_m']-v.chi_B_m)/(v.K_m02_anio*A0**theta)/1000 if usable else np.nan
                    valoreskp.append({'rid':rid,'rio':v.rio,'id_medicion':v.id_medicion,'id_B':v.id_B_corto,'mis':v.mis,
                        'caso_B':case,'id_kp':k['id_kp'],'K_m02_anio':v.K_m02_anio,'chi_B_m':v.chi_B_m,
                        **k,'tau_KP_ka':t,'estado':'utilizable' if usable else 'KP fuera del perfil o aguas abajo de B'})
                    if usable:az.scatter([t],[k['z_kp_m']],marker='*',c=color,s=20,alpha=.65,zorder=4)
            m=np.vstack(matrix);low=m.min(axis=0);high=m.max(axis=0)
            label=f'{b.id_B_corto} · MIS {b.mis} · {len(q)} K'
            at.fill_between(x,low,high,color=color,alpha=.2,label=label)
            at.plot(x,low,color=color,lw=.8);at.plot(x,high,color=color,lw=.8)
            az.plot([],[],color=color,label=label)
            az.scatter([b.T_ka],[b.H_m],marker='o',facecolors='white',edgecolors=color,s=50,zorder=5)
            envolventes.append(pd.DataFrame({'rid':rid,'caso_B':case,'id_B':b.id_B_corto,'chi_m':x,
                'tau_min_ka':low,'tau_max_ka':high,'n_escenarios':len(q)}))
        az.set_xlabel('τ relativo a cada B (ka)');az.set_ylabel('Elevación DEM (m)')
        az.set_title('Perfiles transformados con todos los K')
        at.set_xlabel('χ (m), origen común del perfil');at.set_ylabel('τ relativo a cada B (ka)');at.set_title('Rango de τ por B, sin escoger un K')
        for ax in [az,at]:ax.grid(alpha=.2);ax.set_xlim(left=0);ax.set_ylim(bottom=0)
        az.plot([],[],marker='*',ls='',color='black',label='KP manual')
        az.plot([],[],marker='o',mfc='white',ls='',color='black',label='P_T: τ=T por construcción')
        az.legend(fontsize=8);at.legend(fontsize=8)
        fig.suptitle(f'R{int(rid):02d} · {group.rio.iloc[0]} — τ con K propios\nK constante al extender cada escenario; τ no es una datación tectónica independiente',fontsize=12)
        target=folder/f'R{int(rid):02d}_tau_perfiles_rangos.png';fig.savefig(target,dpi=165,facecolor='white');plt.close(fig);figures.append(target)
    curvas=pd.concat(curvas,ignore_index=True);env=pd.concat(envolventes,ignore_index=True)
    kptable=pd.DataFrame(valoreskp)
    curvas.to_csv(folder/'Tau_todas_las_curvas.csv.gz',index=False,compression='gzip')
    env.to_csv(folder/'Tau_envolventes_por_B.csv.gz',index=False,compression='gzip')
    kptable.to_csv(folder/'Tau_KP_por_escenario.csv',index=False,encoding='utf-8-sig')
    assert max(verificaciones)<1e-8
    return figures,kptable,env,{'max_error_tau_PT_menos_T_ka':max(verificaciones),'n_curvas':resultados.id_medicion.nunique(),'n_nodos_curvas':len(curvas)}
