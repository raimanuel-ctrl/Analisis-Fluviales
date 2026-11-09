
import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
def guardar_png_pdf(fig,path):
    path=pa.version_path(path)
    while path.with_suffix('.pdf').exists():
        stem=path.stem.rsplit('_v',1)[0]
        version=int(path.stem.rsplit('_v',1)[1])+1
        path=path.with_name(f'{stem}_v{version}{path.suffix}')
    fig.savefig(path,dpi=300,bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight')
    return path

import memoria_comun as pa
"""Cálculos y figuras para incorporar en el notebook nuevo, sin ejecutar al importar."""
from pathlib import Path
import json, hashlib, re
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.stats import linregress
from pyproj import Transformer
import matplotlib.pyplot as plt
import plotly.graph_objects as go

COLS_PERFIL = ['distancia_m', 'elevacion_m', 'area_m2', 'x_utm', 'y_utm']

def huella_archivo(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def cargar_perfiles(ruta_manifest, theta, A0):
    if not 0 < theta < 1 or A0 <= 0:
        raise ValueError('Revisar theta y área de referencia.')
    mf = json.loads(Path(ruta_manifest).read_text(encoding='utf-8-sig'))
    if mf.get('estado') != 'completo':
        raise ValueError('La corrida de perfiles no está completa.')
    perfiles, filas = {}, []
    transform = Transformer.from_crs('EPSG:32718', 'EPSG:4326', always_xy=True)
    for rid in range(1, 19):
        item = mf['files'][f'rio_{rid}']
        path = Path(item['path'])
        if not path.is_absolute(): path = Path(ruta_manifest).resolve().parent / path
        if huella_archivo(path) != item['sha256']:
            raise ValueError(f'R{rid:02}: el perfil cambió respecto al manifest.')
        df = pd.read_csv(path, float_precision='round_trip')
        data = df[COLS_PERFIL].to_numpy(float)
        s, z, a = (df[c].to_numpy(float) for c in COLS_PERFIL[:3])
        if not np.isfinite(data).all() or np.any(np.diff(s) <= 0) or np.any(a <= 0):
            raise ValueError(f'R{rid:02}: datos no finitos, distancia no creciente o área inválida.')
        # No se ordena ni se suaviza z: preservar el perfil del DEM original.
        J = cumulative_trapezoid(a ** (-theta), x=s, initial=0.0)
        chi = A0 ** theta * J
        df = df.copy()
        df['chi_cal_m'] = chi
        lon, lat = transform.transform(df.x_utm.iloc[0], df.y_utm.iloc[0])
        perfiles[rid] = dict(df=df, s=s, z=z, a=a, chi=chi, J=J,
                             rio=mf['rios'][str(rid)]['nombre'], latitud=lat,
                             longitud=lon, hash_archivo=item['sha256'], path=str(path),
                             theta=theta,A0=A0)
        filas.append(dict(rid=rid, rio=perfiles[rid]['rio'], latitud=lat,
                          z_inicio_m=z[0], longitud_km=(s[-1]-s[0])/1000,
                          n_nodos=len(s), pasos_z_negativos=int(np.sum(np.diff(z)<0)),
                          sha256_verificado=True))
    return perfiles, pd.DataFrame(filas), mf

def verificar_kp(path, perfiles):
    raw = pd.read_csv(path, sep=';', decimal=',', dtype={'ID': str}).dropna(how='all')
    raw.columns = raw.columns.str.strip()
    necesarios = {'ID', 'id_kp', 'x_utm', 'y_utm', 'elevacion_m', 'perfil_hash'}
    if necesarios - set(raw):
        raise ValueError(f'Faltan columnas KP: {necesarios-set(raw)}')
    filas, vistos = [], set()
    for _, r in raw.iterrows():
        ident = '' if pd.isna(r.id_kp) else str(r.id_kp).strip()
        if pd.isna(r.ID) and ident:
            raise ValueError('KP seleccionado sin ID de río.')
        rid = int(float(r['id'] if pd.isna(r.ID) else r.ID))
        if not ident:
            if r[['x_utm','y_utm','elevacion_m']].notna().any():
                raise ValueError(f'R{rid}: coordenadas sin identificación KP.')
            filas.append(dict(rid=rid, estado='sin_KP_manual', id_kp='', s_kp_m=np.nan,
                               z_kp_m=np.nan, chi_kp_m=np.nan))
            continue
        if ident in vistos or any(v['rid']==rid for v in filas):
            raise ValueError(f'KP duplicado o más de una selección para R{rid}.')
        vistos.add(ident)
        p = perfiles[rid]
        xy = r[['x_utm','y_utm']].to_numpy(float)
        if not np.isfinite(xy).all() or not np.isfinite(float(r.elevacion_m)):
            raise ValueError(f'R{rid}: coordenadas KP no finitas.')
        dist = np.linalg.norm(p['df'][['x_utm','y_utm']].to_numpy(float)-xy, axis=1)
        indices = np.flatnonzero(dist <= 0.10)
        if len(indices) != 1:
            raise ValueError(f'R{rid}: no existe un nodo único a menos de 0.10 m del KP.')
        j = int(indices[0]); dz = abs(p['z'][j]-float(r.elevacion_m))
        if dz > 0.10:
            raise ValueError(f'R{rid}: elevación KP incompatible ({dz:.3f} m).')
        phash = hashlib.sha256(np.asarray(p['df'][COLS_PERFIL].to_numpy(float),
                                         dtype='<f8', order='C').tobytes()).hexdigest()
        filas.append(dict(rid=rid, id_kp=ident, estado='verificado_por_XY_z',
                          s_kp_m=p['s'][j], z_kp_m=p['z'][j], chi_kp_m=p['chi'][j],
                          indice_actual=j, distancia_xy_m=dist[j], diferencia_z_m=dz,
                          hash_original=r.perfil_hash, hash_actual=phash,
                          s_original_m=r.distancia_m))
    return pd.DataFrame(filas)

def inventario_desde_fuentes(old, terrazas, asignaciones, alturas_propias=None):
    """Conserva las asociaciones elegidas, pero lee H, errores y MIS de su fuente."""
    alturas_propias = alturas_propias or {}
    if old.rid.duplicated().any() or terrazas.terrace_id.duplicated().any():
        raise ValueError('Hay identificadores duplicados en las entradas.')
    tf = terrazas.set_index('terrace_id')
    edades_muestra = {'CONT 11': (124.1,10.8), 'PUHE 11': (111.8,9.5), 'PUHE 10': (102.3,9.3)}
    filas=[]
    for rid in range(1,19):
        r=old.loc[old.rid==rid].iloc[0]
        tid=asignaciones.get(rid)
        fila=dict(rid=rid, rio=r.rio, terraza_id=tid or '', H_original_m=r.cota_T_m,
                  edad_original_ka=r.edad_ka, OD_original_pct=r.od_pct,
                  K_anterior_guardado=r.K_ols, K_George_guardado=r.K_george,
                  U_original_guardado=r.U_terraza_m_ka,
                  id_KP_original=r.id_kp, s_KP_original_m=r.s_kp_m,
                  asociacion='sin marcador seleccionado' if tid is None else 'correlación geomorfológica pendiente de revisión local',
                  H_m=np.nan, error_H_publicado_m=np.nan, edad_MIS_ka=np.nan,
                  error_edad_MIS_ka=np.nan, U_JM_m_ka=np.nan, error_U_JM_m_ka=np.nan,
                  MIS='', muestra_referencia='', edad_muestra_ka=np.nan,
                  error_edad_muestra_ka=np.nan, latitud_marcador=np.nan)
        if tid:
            if tid not in tf.index:
                raise ValueError(f'R{rid}: terraza no encontrada: {tid}')
            t=tf.loc[tid]
            muestra='CONT 11' if rid in [3,4,5] else ('PUHE 10' if rid in [10,11] else 'PUHE 11')
            fila.update(H_m=float(t.elevacion_m), error_H_publicado_m=float(t.error_elevacion_m),
                        MIS=t.mis_asignado, edad_MIS_ka=float(t.edad_tabla_2_articulo_ka),
                        error_edad_MIS_ka=float(t.incertidumbre_edad_ka),
                        muestra_referencia=muestra, edad_muestra_ka=edades_muestra[muestra][0],
                        error_edad_muestra_ka=edades_muestra[muestra][1],
                        U_JM_m_ka=float(t.tasa_alzamiento_m_ka),
                        error_U_JM_m_ka=float(t.error_tasa_alzamiento_m_ka),
                        latitud_marcador=float(t.latitud_wgs84),
                        fuente_H=t.archivo_fuente, fila_fuente=int(t.fila_fuente_xls),
                        nota_cronologia=t.nota_edad_referencia,
                        z_mar_fuente_U_JM=float(t.elevacion_m-t.tasa_alzamiento_m_ka*t.edad_referencia_ka),
                        edad_fuente_U_JM_ka=float(t.edad_referencia_ka))
            if rid in alturas_propias:
                H,eH=alturas_propias[rid]
                fila.update(H_m=float(H),error_H_publicado_m=float(eH),fuente_H='altura propia indicada por usuario')
        filas.append(fila)
    inv=pd.DataFrame(filas)
    inv.attrs['firma_selecciones']=json.dumps({'terrazas':asignaciones,'alturas':alturas_propias},sort_keys=True)
    return inv

def cruces_cota(p, H):
    """Intersecciones lineales con todos los segmentos, sin invertir una z no monótona."""
    filas=[]
    s,z,chi,J=p['s'],p['z'],p['chi'],p['J']
    if not np.isfinite(H): return pd.DataFrame(columns=['cruce','s_PT_m','chi_PT_m','J_PT','Mchi_local','direccion','meseta'])
    for j in range(len(s)-1):
        dz=z[j+1]-z[j]
        if dz==0:
            if np.isclose(z[j],H,atol=1e-9,rtol=0):
                filas.append(dict(s_PT_m=s[j],chi_PT_m=chi[j],J_PT=J[j],
                                  Mchi_local=0.,direccion='meseta',meseta=True,indice=j))
            continue
        if not (min(z[j:j+2]) <= H <= max(z[j:j+2])): continue
        f=(H-z[j])/dz
        punto=s[j]+f*(s[j+1]-s[j])
        if filas and np.isclose(punto,filas[-1]['s_PT_m'],atol=1e-6,rtol=0): continue
        filas.append(dict(s_PT_m=punto,chi_PT_m=chi[j]+f*(chi[j+1]-chi[j]),
                          J_PT=J[j]+f*(J[j+1]-J[j]),
                          Mchi_local=dz/(chi[j+1]-chi[j]),
                          direccion='ascendente' if dz>0 else 'descendente',meseta=False,indice=j))
    df=pd.DataFrame(filas)
    if not df.empty: df.insert(0,'cruce',np.arange(1,len(df)+1))
    return df

def ajuste_hasta(p, s_fin, min_nodos=5):
    mask=(p['s']>=p['s'][0]) & (p['s']<=s_fin)
    x,z=p['chi'][mask],p['z'][mask]
    if len(x)<min_nodos or np.ptp(x)<=0:
        return dict(Mchi=np.nan,ksn=np.nan,R2=np.nan,RMSE_m=np.nan,
                    intercepto_m=np.nan,n_fit=len(x),sigma_Mchi_formal=np.nan)
    reg=linregress(x,z)
    res=z-(reg.intercept+reg.slope*x)
    return dict(Mchi=reg.slope,R2=reg.rvalue**2,RMSE_m=np.sqrt(np.mean(res**2)),
                intercepto_m=reg.intercept,n_fit=len(x),sigma_Mchi_formal=reg.stderr)

def calcular_tabla(inventario, perfiles, kp, cfg):
    A0,theta=cfg['A0'],cfg['theta']; factor=A0**theta
    if cfg['n']!=1: raise ValueError('Estas ecuaciones requieren n=1.')
    filas, todos=[] , []
    for _,r in inventario.iterrows():
        rid=int(r.rid);p=perfiles[rid]
        if p['theta']!=theta or p['A0']!=A0:
            raise ValueError('Cambió theta/A0: recargar perfiles y verificar KP antes de calcular.')
        o=r.to_dict()
        o.update(latitud=p['latitud'],z_inicio_m=p['z'][0],chi_B_m=float(cfg['chi_B_por_rio'].get(rid,0.)),
                 error_chi_B_m=float(cfg['error_chi_B_por_rio'].get(rid,0.)),
                 K_George=np.nan,error_K_George=np.nan,K_pendiente=np.nan,error_K_pendiente=np.nan,
                 K_JM_pendiente=np.nan,error_K_JM_pendiente=np.nan,
                 U_terraza=np.nan,error_U_terraza=np.nan,U_George=np.nan,error_U_George=np.nan,
                 s_PT_m=np.nan,chi_PT_m=np.nan,delta_chi_m=np.nan,Mchi=np.nan,ksn=np.nan,
                 R2=np.nan,RMSE_m=np.nan,intercepto_m=np.nan,n_fit=0,t_viaje_KP_ka=np.nan,
                 n_cruces=0,cruce_elegido=np.nan,estado_K='sin marcador seleccionado',
                 estado_error_K='no calculable',error_DEM_m=0.,nodos_nivel_no_monotono=True)
        o['estado_K_pendiente']='sin tramo de calibración'
        sub=kp[kp.rid==rid]
        if len(sub):
            o.update(id_KP=sub.iloc[0].id_kp,estado_KP=sub.iloc[0].estado,
                     s_KP_m=sub.iloc[0].s_kp_m,chi_KP_m=sub.iloc[0].chi_kp_m)
        else: o.update(id_KP='',estado_KP='sin_KP_manual',s_KP_m=np.nan,chi_KP_m=np.nan)
        if not np.isfinite(r.H_m): filas.append(o); continue
        if cfg['edad_modo']=='MIS':
            T,eT=r.edad_MIS_ka,r.error_edad_MIS_ka
        elif cfg['edad_modo']=='muestra':
            T,eT=r.edad_muestra_ka,r.error_edad_muestra_ka
        else: raise ValueError('edad_modo debe ser MIS o muestra.')
        sea=cfg['nivel_mar_por_MIS'][r.MIS]
        SL,eSL=float(sea['valor_m']),float(sea['error_publicado_m'])
        eH=float(r.error_H_publicado_m)
        if not np.isfinite(T) or T<=0 or not np.isfinite(eT) or eT<0:
            raise ValueError(f'R{rid}: edad o error de edad inválidos.')
        if not np.isfinite(eH) or eH<0: raise ValueError(f'R{rid}: falta error de altura.')
        # Se propagan magnitudes publicadas sin afirmar un IC uniforme entre fuentes.
        U=(r.H_m-SL)/T
        eU=np.sqrt((eH/T)**2+(eSL/T)**2+((r.H_m-SL)*eT/T**2)**2)
        o.update(edad_ka=T,error_edad_ka=eT,z_mar_m=SL,error_z_mar_m=eSL,
                 U_terraza=U,error_U_terraza=eU)
        cr=cruces_cota(p,r.H_m)
        o['n_cruces']=len(cr)
        if not cr.empty:
            cr=cr.assign(rid=rid,H_m=r.H_m,K_candidato=(cr.chi_PT_m-o['chi_B_m'])/(factor*T*1000))
            todos.append(cr)
        if cr.empty:
            o['estado_K']='cota bajo perfil disponible' if r.H_m<np.min(p['z']) else 'cota fuera de rango / sin cruce'
            filas.append(o); continue
        elegido=cfg['cruce_por_rio'].get(rid)
        if elegido is None:
            ascendentes=cr.loc[(cr.direccion=='ascendente') & ~cr.meseta]
            if ascendentes.empty:
                o['estado_K']='sin cruce ascendente; revisar';filas.append(o);continue
            seleccionado=ascendentes.iloc[0]
        else:
            candidatos=cr.loc[cr.cruce==int(elegido)]
            if len(candidatos)!=1: raise ValueError(f'R{rid}: cruce solicitado no existe.')
            seleccionado=candidatos.iloc[0]
        c=seleccionado
        delta=float(c.chi_PT_m-o['chi_B_m'])
        o.update(s_PT_m=c.s_PT_m,chi_PT_m=c.chi_PT_m,delta_chi_m=delta,
                 cruce_elegido=int(c.cruce),Mchi_local=c.Mchi_local)
        if delta<=0 or bool(c.meseta):
            o['estado_K']='viaje no positivo o cruce sobre meseta';filas.append(o);continue
        KG=delta/(factor*T*1000)
        o.update(K_George=KG, c_chi_m_ka=delta/T,
                 estado_K='calculado condicional; cruce único' if len(cr)==1 else 'calculado condicional; cruce múltiple')
        # A0 cancela con el mismo esquema de interpolación. B!=0 incluye J_B=chi_B/A0^m.
        KI=(c.J_PT-o['chi_B_m']/factor)/(T*1000)
        o['K_integral']=KI
        o['error_relativo_identidad']=abs(KG-KI)/max(abs(KG),1e-30)
        Mloc=abs(float(c.Mchi_local))
        if Mloc>0:
            eDelta=np.sqrt((eH/Mloc)**2+o['error_chi_B_m']**2)
            eKG=np.sqrt((eDelta/(factor*T*1000))**2+(KG*eT/T)**2)
            o.update(error_K_George=eKG,error_delta_chi_m=eDelta,
                     estado_error_K='linealizado en cruce elegido; B y DEM fijos' if o['error_chi_B_m']==0 else 'linealizado; incertidumbre B indicada')
            # Rango de cruces alternativos y errores analíticos quedan separados.
            o['aproximacion_error_fragil']=(eH>=abs(r.H_m-p['z'][0]) or eKG>=KG or len(cr)>1)
        fit=ajuste_hasta(p,c.s_PT_m,cfg['min_nodos_fit']);o.update(fit)
        M=fit['Mchi'];o['ksn']=M*factor if np.isfinite(M) else np.nan
        if np.isfinite(M) and M>0:
            o.update(K_pendiente=(U/1000)/(M*factor),error_K_pendiente=(eU/1000)/(M*factor),
                     U_George=KG*factor*M*1000,
                     error_U_George=o['error_K_George']*factor*M*1000,
                     K_JM_pendiente=(r.U_JM_m_ka/1000)/(M*factor),
                     error_K_JM_pendiente=(r.error_U_JM_m_ka/1000)/(M*factor))
            o['estado_K_pendiente']='calculado condicional; ajuste considerado fijo'
        else:
            o['estado_K_pendiente']='pendiente OLS no positiva' if np.isfinite(M) else 'nodos insuficientes para OLS'
        if np.isfinite(o['chi_KP_m']) and o['chi_KP_m']>o['chi_B_m']:
            o['t_viaje_KP_ka']=(o['chi_KP_m']-o['chi_B_m'])/(KG*factor)/1000
        filas.append(o)
    resultados=pd.DataFrame(filas)
    resultados.attrs['cfg_calculo']=json.dumps(cfg,sort_keys=True)
    cruces=pd.concat(todos,ignore_index=True) if todos else pd.DataFrame(columns=['rid','cruce','s_PT_m','chi_PT_m','K_candidato'])
    return resultados,cruces

def comprobar_cfg(df,cfg):
    if df.attrs.get('cfg_calculo')!=json.dumps(cfg,sort_keys=True):
        raise ValueError('Cambió la configuración: ejecutar de nuevo la celda de cálculo antes de graficar/exportar.')

def comparar_original(df,old,perfiles):
    comparacion=df.merge(old[['rid','K_george','K_ols','s_PT_m','chi_PT_m','c_chi_m_yr']],on='rid',suffixes=('','_original'))
    comparacion['K_desde_formula_original']=comparacion.c_chi_m_yr/(1e6**.4)
    comparacion['K_guardado_no_coincide_formula_original']=~np.isclose(comparacion.K_george,comparacion.K_desde_formula_original,rtol=1e-7,atol=1e-14,equal_nan=True)
    comparacion['KP_distancia_cambio_m']=comparacion.s_KP_m-comparacion.s_KP_original_m
    comparacion['razon_K_nuevo_original']=comparacion.K_George/comparacion.K_George_guardado
    return comparacion

ESTILOS_SERIES={
    'George':dict(color='#235789',marker='o',dash='solid'),
    'propio':dict(color='#d67916',marker='s',dash='dash'),
    'JM':dict(color='#7855a1',marker='D',dash='dot'),
    'anterior':dict(color='#777777',marker='x',dash='dashdot'),
}

def datos_series(df,variable):
    if variable=='U':
        return [('George','U_George','error_U_George','U fluvial condicionado a K George'),
                ('propio','U_terraza','error_U_terraza','U de terraza (H − mar) / T'),
                ('JM','U_JM_m_ka','error_U_JM_m_ka','U de J-M en marcador asociado'),
                ('anterior','U_original_guardado',None,'U anterior guardado')]
    return [('George','K_George','error_K_George','K George: viaje en χ'),
            ('propio','K_pendiente','error_K_pendiente','K por pendiente: U terraza / ksn'),
            ('JM','K_JM_pendiente','error_K_JM_pendiente','K por pendiente con U de J-M'),
            ('anterior','K_anterior_guardado',None,'K anterior guardado')]



def perfiles_calibracion(df,perfiles,cruces,rios,carpeta=None):
    figuras=[]
    for rid in rios:
        row=df.loc[df.rid==int(rid)].iloc[0];p=perfiles[int(rid)]
        fig,grid=plt.subplots(2,2,figsize=(12.5,7.8),layout='constrained')
        axes=[grid[0,0],grid[1,0],grid[1,1]]
        zoom=grid[0,1]
        H=row.H_m
        axes[0].plot(p['s']/1000,p['z'],color='#33465a',lw=.9)
        axes[1].plot(p['chi'],p['z'],color='#33465a',lw=.9)
        zoom.plot(p['s']/1000,p['z'],color='#33465a',lw=.9)
        if np.isfinite(H):
            for ax in axes[:2]+[zoom]:ax.axhline(H,color='#aa752b',ls='--',lw=1,label=f'H={H:.2f} m')
            candidates=cruces[cruces.rid==int(rid)]
            for _,c in candidates.iterrows():
                for ax,x in [(axes[0],c.s_PT_m/1000),(zoom,c.s_PT_m/1000),(axes[1],c.chi_PT_m)]:
                    ax.scatter(x,H,s=24,facecolor='white',edgecolor='#aa752b',zorder=4)
            if np.isfinite(row.K_George):
                axes[1].scatter(row.chi_PT_m,H,color='#235789',s=40,zorder=5,label='P_T elegido')
                axes[0].axvspan(p['s'][0]/1000,row.s_PT_m/1000,color='#235789',alpha=.10,label='Tramo del ajuste')
                zoom.axvspan(p['s'][0]/1000,row.s_PT_m/1000,color='#235789',alpha=.10)
                mask=p['s']<=row.s_PT_m
                if np.isfinite(row.Mchi):
                    x=p['chi'][mask];res=p['z'][mask]-(row.intercepto_m+row.Mchi*x)
                    axes[1].plot(x,row.intercepto_m+row.Mchi*x,color='#d67916',lw=1.5,label='Recta OLS libre')
                    axes[2].plot(x,res,color='#33465a',lw=.8);axes[2].axhline(0,color='#888888',lw=.6)
        if np.isfinite(row.s_KP_m):
            zkp=np.interp(row.s_KP_m,p['s'],p['z'])
            axes[0].scatter(row.s_KP_m/1000,zkp,marker='*',color='#8b2d3d',s=90,label='KP manual')
            axes[1].scatter(row.chi_KP_m,zkp,marker='*',color='#8b2d3d',s=90)
        axes[0].set(xlabel='Distancia desde outlet (km)',ylabel='Elevación (m)',title='Perfil completo y cota marina')
        axes[1].set(xlabel='χ (m)',ylabel='Elevación (m)',title='Calibración en χ')
        zoom.set(xlabel='Distancia desde outlet (km)',ylabel='Elevación (m)',title='Zoom del tramo inferior y sus cruces')
        if np.isfinite(row.s_PT_m):
            # Mostrar los cruces del tramo costero, incluso las alternativas posteriores.
            csub=cruces[cruces.rid==int(rid)]
            fin=min(p['s'][-1], max(float(row.s_PT_m)*1.3,
                                    float(csub.s_PT_m.max())*1.1 if len(csub) else 0,
                                    p['s'][0]+250))
            near=p['s']<=fin
            zz=p['z'][near]
            margin=max(3.,(max(np.max(zz),H)-min(np.min(zz),H))*.10)
            zoom.set_xlim(p['s'][0]/1000,fin/1000)
            zoom.set_ylim(min(np.min(zz),H)-margin,max(np.max(zz),H)+margin)
            axes[1].set_xlim(p['chi'][0],float(np.interp(fin,p['s'],p['chi'])))
            axes[1].set_ylim(zoom.get_ylim())
        axes[2].set(xlabel='χ (m)',ylabel='Residual (m)',title='¿El tramo se aproxima a una recta?')
        for ax in axes+[zoom]:ax.grid(alpha=.15)
        axes[0].legend(fontsize=7,loc='best');axes[1].legend(fontsize=7,loc='best')
        fig.suptitle(f"R{int(rid):02} · {row.rio} · {row.estado_K}",fontsize=12)
        if carpeta:guardar_png_pdf(fig,Path(carpeta)/f'{int(rid):02}_Calibracion.png')
        figuras.append(fig)
    return figuras

def diagnostico_cruces(df,cruces,carpeta=None):
    fig,axes=plt.subplots(1,2,figsize=(12.8,4.6),layout='constrained')
    axes[0].axhline(0,color='#999999',lw=.8)
    v=df[np.isfinite(df.H_m)].sort_values('latitud',ascending=False)
    axes[0].bar(v.rid.astype(str),v.H_m-v.z_inicio_m,color=np.where(v.H_m-v.z_inicio_m>0,'#235789','#a63945'))
    axes[0].set(xlabel='Río (ID)',ylabel='H terraza − z inicio del perfil (m)',title='Disponibilidad del tramo inferior')
    for rid,g in cruces.groupby('rid'):
        y=g.K_candidato.to_numpy(float); good=y>0
        axes[1].scatter(np.repeat(rid,np.sum(good)),y[good],facecolor='white',edgecolor='#aa752b',s=38)
    good=df.K_George.notna()
    axes[1].scatter(df.loc[good,'rid'],df.loc[good,'K_George'],color='#235789',s=25,label='Cruce elegido',zorder=5)
    axes[1].set(xlabel='Río (ID)',ylabel='K candidato (m^(1−2θ)/año)',title='Cruces alternativos: sensibilidad geométrica',yscale='log')
    axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.15)
    if carpeta:guardar_png_pdf(fig,Path(carpeta)/'diagnostico_cruces.png')
    return fig

def sensibilidad_B(df,theta,A0,rios,rango=(-200,200),carpeta=None):
    fig,ax=plt.subplots(figsize=(9.5,4.5),layout='constrained')
    bs=np.linspace(*rango,151)
    for rid in rios:
        r=df.loc[df.rid==int(rid)].iloc[0]
        if not np.isfinite(r.K_George):continue
        K=(r.chi_PT_m-bs)/(A0**theta*r.edad_ka*1000)
        K[K<=0]=np.nan
        ax.plot(bs,K,label=f'R{int(rid):02} · {r.rio}')
    ax.axvline(0,color='#777777',ls='--',lw=.8)
    ax.set(xlabel='χ(B) propuesto (m de χ)',ylabel='K (m^(1−2θ)/año)',yscale='log',
           title='Sensibilidad al origen B · rango exploratorio, sin reconstrucción geométrica')
    ax.legend(fontsize=8);ax.grid(alpha=.15)
    if carpeta:guardar_png_pdf(fig,Path(carpeta)/'sensibilidad_origen_B.png')
    return fig


"""Extensión explícita de estimaciones regionales y comparación con Grant."""
from matplotlib.lines import Line2D

def agregar_K_regional(resultados, perfiles, kp, terrazas, cfg, cfg_regional):
    comprobar_cfg(resultados,cfg)
    referencias=terrazas.loc[terrazas.mis_asignado.eq(cfg_regional['MIS'])].copy()
    for col in ['latitud_wgs84','tasa_alzamiento_m_ka','error_tasa_alzamiento_m_ka']:
        referencias=referencias.loc[np.isfinite(pd.to_numeric(referencias[col],errors='coerce'))]
    referencias=referencias.sort_values('terrace_id')
    if referencias.empty:
        raise ValueError('No hay referencias de uplift regional compatibles con el MIS elegido.')
    tabla=resultados.copy()
    factor=cfg['A0']**cfg['theta']
    iniciales=set(cfg_regional['rios_iniciales'])
    filas=[]
    for _,r in tabla.iterrows():
        rid=int(r.rid);p=perfiles[rid]
        referencia_elegida=cfg_regional['referencia_por_rio'].get(rid)
        if referencia_elegida:
            candidatas=referencias.loc[referencias.terrace_id.eq(referencia_elegida)]
            if len(candidatas)!=1:
                raise ValueError(f'R{rid}: referencia regional ausente o incompatible con MIS.')
            t=candidatas.iloc[0];tipo='referencia regional indicada'
        elif cfg_regional.get('conservar_asociadas',True) and r.terraza_id in set(referencias.terrace_id):
            t=referencias.loc[referencias.terrace_id.eq(r.terraza_id)].iloc[0]
            tipo='marcador del MIS regional previamente seleccionado'
        else:
            dist_lat=np.abs(referencias.latitud_wgs84.to_numpy(float)-p['latitud'])
            t=referencias.iloc[int(np.argmin(dist_lat))]
            tipo='marcador MIS 5e más próximo en latitud' if cfg_regional['MIS']=='MIS 5e' else 'marcador del MIS elegido más próximo en latitud'
        dxy=np.hypot(float(t.este_utm18s_m)-p['df'].x_utm.iloc[0],
                     float(t.norte_utm18s_m)-p['df'].y_utm.iloc[0])/1000.
        seleccion=kp.loc[(kp.rid==rid)&kp.s_kp_m.notna()]
        if len(seleccion)==1:
            limite=float(seleccion.iloc[0].s_kp_m)
            # Excluir el propio nodo KP, de acuerdo con el ajuste bajo KP.
            fit=ajuste_hasta(p,np.nextafter(limite,-np.inf),cfg['min_nodos_fit'])
            tramo='outlet–KP manual; nodo KP excluido'
        elif rid==7 and np.isfinite(r.s_PT_m):
            limite=float(r.s_PT_m)
            fit=ajuste_hasta(p,limite,cfg['min_nodos_fit'])
            tramo='outlet–P_T; Chanco sin KP manual'
        else:
            raise ValueError(f'R{rid}: falta un tramo definido para la estimación regional.')
        M=fit['Mchi'];ksn=M*factor
        if not np.isfinite(ksn) or ksn<=0:
            raise ValueError(f'R{rid}: pendiente no positiva en el tramo regional seleccionado.')
        U=float(t.tasa_alzamiento_m_ka);eU=float(t.error_tasa_alzamiento_m_ka)
        if U<=0 or eU<0:
            raise ValueError(f'R{rid}: la referencia regional seleccionada no tiene U positivo o error válido.')
        filas.append(dict(rid=rid,grupo_rios='seis iniciales' if rid in iniciales else 'ríos adicionales',
            K_regional=(U/1000.)/ksn,error_K_regional=(eU/1000.)/ksn,
            U_regional=U,error_U_regional=eU,ksn_regional=ksn,Mchi_regional=M,
            R2_regional=fit['R2'],RMSE_regional_m=fit['RMSE_m'],n_fit_regional=fit['n_fit'],
            limite_regional_m=limite,tramo_regional=tramo,
            referencia_U_regional=t.terrace_id,latitud_U_regional=float(t.latitud_wgs84),
            delta_latitud_ref_grados=abs(float(t.latitud_wgs84)-p['latitud']),
            distancia_referencia_outlet_km=dxy,MIS_U_regional=t.mis_asignado,
            metodo_U_regional=tipo,estado_regional='estimación condicional con U transferido; DEM y ajuste fijos'))
    extra=pd.DataFrame(filas).set_index('rid')
    for col in extra.columns:tabla[col]=tabla.rid.map(extra[col])
    tabla.attrs=resultados.attrs.copy()
    tabla.attrs['cfg_regional']=json.dumps(cfg_regional,sort_keys=True)
    return tabla

def validar_regional(tabla,cfg_regional):
    if tabla.attrs.get('cfg_regional')!=json.dumps(cfg_regional,sort_keys=True):
        raise ValueError('Cambió la configuración regional: recalcular la extensión antes de graficar.')

def guardar_regional(tabla,cfg_regional,carpeta):
    validar_regional(tabla,cfg_regional)
    columnas=['rid','rio','grupo_rios','latitud','K_regional','error_K_regional','ksn_regional',
              'U_regional','error_U_regional','referencia_U_regional','latitud_U_regional',
              'distancia_referencia_outlet_km','delta_latitud_ref_grados','tramo_regional',
              'limite_regional_m','n_fit_regional','R2_regional','RMSE_regional_m','estado_regional']
    tabla[columnas].to_csv(Path(carpeta)/'K_regional_18_rios.csv',index=False,encoding='utf-8-sig')
    (Path(carpeta)/'configuracion_regional.json').write_text(json.dumps(cfg_regional,indent=2,ensure_ascii=False),encoding='utf-8')

_datos_series_anteriores=datos_series
def datos_series(df,variable):
    series=_datos_series_anteriores(df,variable)
    if variable=='K' and 'K_regional' in df:
        series.append(('regional','K_regional','error_K_regional','K regional: U J-M / ksn bajo KP'))
    return series

ESTILOS_SERIES['regional']=dict(color='#956244',marker='^',dash='dashdot')

def validar_franja_grant(variable,opciones,cfg):
    if variable=='K' and opciones.get('Grant',False):
        if cfg is None or not np.isclose(cfg['theta'],.40) or cfg['n']!=1:
            raise ValueError('Los K de Grant de esta franja corresponden exclusivamente a theta=0.40, n=1.')
        if not (0<opciones['rango_Grant'][0]<=opciones['rango_Grant'][1]):
            raise ValueError('Rango de Grant inválido.')

def grupos_simbolos(v,opciones):
    ids=set(opciones.get('rios_iniciales',[3,4,7,8,9,10]))
    return v.rid.isin(ids).to_numpy(bool)

def grafico_latitud(df,variable,opciones,carpeta=None,cfg=None):
    if cfg is not None:comprobar_cfg(df,cfg)
    if variable=='K' and opciones['series'].get('regional',{}).get('mostrar',False):
        validar_regional(df,CFG_REGIONAL)
    validar_franja_grant(variable,opciones,cfg)
    v=df.sort_values('latitud',ascending=False)
    inicial=grupos_simbolos(v,opciones)
    fig,ax=plt.subplots(figsize=(13.8,6.5),layout='constrained')
    if opciones.get('Grant',False):
        lo,hi=opciones['rango_Grant']
        label=opciones.get('etiqueta_Grant',f'Grant: {lo:.4g}–{hi:.4g}')
        ax.axhspan(lo,hi,color='#76b5a2',alpha=.22,label=label,zorder=0)
        ax.axhline(lo,color='#428a72',lw=.8,ls=':');ax.axhline(hi,color='#428a72',lw=.8,ls=':')
    for clave,col,err,label in datos_series(v,variable):
        opt=opciones['series'].get(clave,{'mostrar':False})
        if not opt.get('mostrar',False):continue
        style=ESTILOS_SERIES[clave];y=v[col].to_numpy(float)
        if opt.get('linea',True):
            ax.plot(v.latitud,y,color=style['color'],linestyle={'solid':'-','dash':'--','dot':':','dashdot':'-.'}[style['dash']],lw=1.3,alpha=.70,zorder=2)
        if opt.get('barras',True) and err:
            mask=np.isfinite(y)&np.isfinite(v[err].to_numpy(float))
            ax.errorbar(v.loc[mask,'latitud'],y[mask],yerr=v.loc[mask,err],fmt='none',ecolor=style['color'],capsize=2.5,lw=.9,alpha=.75,zorder=3)
        good=np.isfinite(y)
        for mask,symbol,abierto in [(good&inicial,opciones.get('simbolo_inicial','o'),False),
                                   (good&~inicial,opciones.get('simbolo_adicional','^'),True)]:
            ax.scatter(v.loc[mask,'latitud'],y[mask],marker=symbol,s=48,
                       facecolors='none' if abierto else style['color'],
                       edgecolors=style['color'],linewidths=1.3,zorder=5)
        ax.plot([],[],color=style['color'],lw=1.4,label=label)
        if opciones.get('etiquetas',True) and clave==opciones.get('serie_etiquetas','George'):
            for x,yy,rid in zip(v.loc[good,'latitud'],y[good],v.loc[good,'rid']):
                ax.annotate(f'R{int(rid):02}',(x,yy),xytext=(3,8),textcoords='offset points',fontsize=8)
    if opciones.get('fallas',False):
        for nombre,lat in opciones.get('latitudes_fallas',{}).items():
            ax.axvline(lat,color='#a42e38',ls='--',lw=1.25,zorder=1)
            ax.text(lat,.97,nombre,transform=ax.get_xaxis_transform(),rotation=90,va='top',fontsize=8,color='#a42e38')
    if variable=='K' and opciones.get('log',True):ax.set_yscale('log')
    ax.set_xlim(df.latitud.max()+.035,df.latitud.min()-.035)
    ax.set_xlabel('Latitud (°, sur negativo)   Norte → Sur')
    ax.set_ylabel('U (mm/año = m/ka)' if variable=='U' else r'K (m$^{1-2\theta}$/año)')
    ax.set_title('Uplift y latitud' if variable=='U' else 'K de los 18 ríos y comparación con Grant',loc='left',fontsize=14)
    ax.grid(alpha=.16)
    handles,labels=ax.get_legend_handles_labels()
    handles.extend([Line2D([],[],marker='o',color='#777777',lw=0,markersize=6),
                    Line2D([],[],marker='^',color='#777777',markerfacecolor='none',lw=0,markersize=6)])
    labels.extend(['Seis ríos iniciales (círculos)','Ríos adicionales (triángulos abiertos)'])
    ax.legend(handles,labels,loc='best',fontsize=8)
    nota='K regional = U MIS 5e transferido / ksn del tramo indicado. DEM y ajuste fijos; referencias en la tabla.' if variable=='K' else 'B=0; errores publicados; U fluvial condicionado a K George.'
    if variable=='K':nota+=' La franja Grant abarca cuatro valores centrales; no es un intervalo de confianza.'
    if opciones.get('fallas',False):nota+='\nLíneas rojas: extremos occidentales de las trazas CHAF.'
    ax.text(0,-.16,nota,transform=ax.transAxes,fontsize=8,color='#555555')
    if carpeta:guardar_png_pdf(fig,Path(carpeta)/f'{variable}_vs_latitud.png')
    return fig

def interactivo_latitud(df,variable,opciones,carpeta,cfg=None):
    if cfg is not None:comprobar_cfg(df,cfg)
    if variable=='K' and opciones['series'].get('regional',{}).get('mostrar',False):
        validar_regional(df,CFG_REGIONAL)
    validar_franja_grant(variable,opciones,cfg)
    v=df.sort_values('latitud',ascending=False);inicial=grupos_simbolos(v,opciones)
    simbolos=np.where(inicial,'circle','triangle-up-open').tolist()
    fig=go.Figure()
    for clave,col,err,label in datos_series(v,variable):
        opt=opciones['series'].get(clave,{'mostrar':False});style=ESTILOS_SERIES[clave]
        custom=[]
        for _,r in v.iterrows():
            texto=f"R{int(r.rid):02} · {r.rio}<br>{r.grupo_rios}<br>"
            if clave=='regional':
                texto+=(f"U transferido: {r.U_regional:.4g} m/ka · {r.MIS_U_regional}<br>"
                        f"Referencia: {r.referencia_U_regional} · lat {r.latitud_U_regional:.4f}°<br>"
                        f"Distancia al outlet: {r.distancia_referencia_outlet_km:.2f} km<br>"
                        f"Tramo: {r.tramo_regional}<br>ksn: {r.ksn_regional:.4g}<br>"
                        f"R²: {r.R2_regional:.3f} · RMSE: {r.RMSE_regional_m:.2f} m<br>Estimación regional condicional")
            else:
                texto+=(f"Terraza: {r.terraza_id or 'sin selección'}<br>{r.estado_K}<br>"
                        f"K por pendiente: {r.estado_K_pendiente}<br>Cruces: {r.n_cruces}")
            if err:texto+=f'<br>Error propagado: {r[err]:.4g}'
            custom.append(texto)
        fig.add_trace(go.Scatter(x=v.latitud,y=v[col],customdata=custom,
            mode='lines+markers' if opt.get('linea',True) else 'markers',name=label,
            legendgroup=clave,visible=True if opt.get('mostrar',False) else 'legendonly',
            connectgaps=False,line=dict(color=style['color'],dash=style['dash'],width=1.6),
            marker=dict(size=9,color=style['color'],symbol=simbolos,line=dict(width=1.4)),
            error_y=dict(type='data',array=v[err] if err else np.full(len(v),np.nan),
                         visible=bool(err and opt.get('barras',True)),thickness=1,width=3),
            hovertemplate='%{customdata}<br>Latitud del río: %{x:.4f}°<br>Valor: %{y:.4g}<extra>%{fullData.name}</extra>'))
    # Franja como traza: se puede apagar desde la leyenda además de True/False.
    if variable=='K':
        lo,hi=opciones['rango_Grant'];xmin=df.latitud.min()-.035;xmax=df.latitud.max()+.035
        valores='<br>'.join(f'{k}: {value:.4e}' for k,value in opciones['valores_Grant'].items())
        fig.add_trace(go.Scatter(x=[xmin,xmax,xmax,xmin,xmin],y=[lo,lo,hi,hi,lo],
            fill='toself',fillcolor='rgba(118,181,162,0.22)',line=dict(width=0),
            name=opciones['etiqueta_Grant'],legendgroup='Grant',
            visible=True if opciones.get('Grant',False) else 'legendonly',
            hovertemplate='Grant · θ=0.40<br>'+valores+'<br>Rango de valores centrales, sin ampliar por error<extra></extra>'))
        fig.data=tuple(fig.data[-1:])+tuple(fig.data[:-1])  # Franja detrás de los puntos.
    for symbol,label in [('circle','Seis ríos iniciales'),('triangle-up-open','Ríos adicionales')]:
        fig.add_trace(go.Scatter(x=[None],y=[None],mode='markers',marker=dict(symbol=symbol,color='#777777',size=8),
                                name=label,hoverinfo='skip',showlegend=True))
    if opciones.get('fallas',False):
        for nombre,lat in opciones.get('latitudes_fallas',{}).items():
            fig.add_vline(x=lat,line_dash='dash',line_color='#a42e38',line_width=1.5,annotation_text=nombre,annotation_textangle=-90,annotation_position='top left')
    indices_series=[i for i,t in enumerate(fig.data) if t.legendgroup in ESTILOS_SERIES]
    fig.update_layout(template='plotly_white',height=690,
        title=('Uplift vs. latitud' if variable=='U' else f"K de los 18 ríos · comparación con Grant · θ={cfg['theta']:.2f}"),
        xaxis=dict(title='Latitud (°, sur negativo) · Norte → Sur',range=[df.latitud.max()+.035,df.latitud.min()-.035]),
        yaxis=dict(title='U (mm/año)' if variable=='U' else 'K (m^(1−2θ)/año)',type='log' if variable=='K' and opciones.get('log',True) else 'linear'),
        legend=dict(orientation='h',y=-.20),margin=dict(t=110,b=170),
        updatemenus=[dict(type='buttons',direction='right',x=0,y=1.16,
            buttons=[dict(label='Mostrar barras',method='restyle',args=[{'error_y.visible':True},indices_series]),
                     dict(label='Ocultar barras',method='restyle',args=[{'error_y.visible':False},indices_series])])])
    path=pa.version_path(Path(carpeta)/f'{variable}_vs_latitud_interactivo.html')
    fig.write_html(path,include_plotlyjs=True,config={'scrollZoom':True,'displaylogo':False})
    return fig,path
