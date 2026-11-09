"""Definición vigente de K, τ y tasas aparentes. No adopta un K único por río."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import json,shutil
import numpy as np
import pandas as pd
import geopandas as gpd
from pyproj import Transformer
from shapely.geometry import Point
import matplotlib.pyplot as plt
import revision_b_base as rb
import escenarios_k_manuales as km
import k_dem_terraza_v8 as v8
import red_fallas as rf
import calibracion_tau as ct
import cronologia_jm2015 as cronologia

LL=Transformer.from_crs(32718,4326,always_xy=True)

def geometria(base,apoyo,out,reconstruir=False,csv_nuevo=None):
    base,apoyo,out=map(Path,[base,apoyo,out]);F=apoyo/'Fuentes'
    anterior=base/'Productos/K_manuales_v7/20261016T223006_387828Z'
    old=pd.read_csv(base/'Productos/03_K_U/v2/resultados_K_U_v1.csv')
    nombres=old.set_index('rid').rio.to_dict()
    if not reconstruir:
        frozen=json.loads((base/'Productos/K_manuales_v8/20261016T224600_778141Z/Trazabilidad_y_verificacion.json').read_text(encoding='utf8'))
        for file,expected in frozen['fuentes_sha256'].items():
            path=Path(file)
            if path.suffix.lower() in ['.tif','.csv'] and rf.sha(path)!=expected:
                raise ValueError('Cambió una fuente congelada de la calibración: '+str(path))
        if csv_nuevo and Path(csv_nuevo).exists() and rf.sha(csv_nuevo)!=rf.sha(F/'Puntos_B_todos_los_rios.csv'):
            raise ValueError('El CSV manual cambió. Activa RECONSTRUIR_B_DESDE_CSV para incorporar los nuevos puntos; no se usarán B antiguos silenciosamente.')
        casos=pd.read_csv(anterior/'B_todos.csv',float_precision='round_trip')
        perfiles={r:pd.read_csv(anterior/f'perfil_R{r:02d}_utilizado.csv',float_precision='round_trip') for r in range(1,19)}
        qc=pd.DataFrame();origen='geometría congelada y verificada; recalcular H, PT y K'
    else:
        ruta=Path(csv_nuevo) if csv_nuevo else F/'Puntos_B_todos_los_rios.csv'
        puntos=km.leer_puntos(ruta)
        perfiles,nombres,_,cfg,_=rb.cargar_perfiles(base,list(range(1,19)))
        for p in perfiles.values():p['tramo_recuperado_aguas_abajo']=False
        for rid,g in puntos.groupby('rid'):
            p=perfiles[int(rid)];line=rb.linea_rio(p)
            if any(line.project(Point(v.x_utm18s_m,v.y_utm18s_m))<1e-6 for _,v in g.iterrows()):
                perfiles[int(rid)],_=km.extender_perfil(base,p,int(rid),8000.)
        cfg=json.loads((base/'Configuracion.json').read_text(encoding='utf-8-sig'))
        casos,qc=km.construir_casos(out,perfiles,puntos,F/'AnalisisK.shp',F/'Exclusiones_otro_experimento.json',base/cfg['dem'],100.)
        shutil.copy2(ruta,out/'Puntos_B_entrada.csv');origen='B reconstruidos desde CSV manual y marcas originales'
    casos.to_csv(out/'B_utilizados.csv',index=False,encoding='utf-8-sig')
    for rid,p in perfiles.items():p.to_csv(out/f'perfil_R{rid:02d}_utilizado.csv',index=False,encoding='utf-8-sig')
    return casos,perfiles,nombres,qc,origen

def calcular_K(casos,perfiles,nombres,catalogo,dem):
    refs=cronologia.aplicar(v8.referencias_dem(casos,catalogo,nombres,dem))
    res,qc=km.calcular(refs,perfiles)
    d=res.sort_values(['rid','caso_B','s_PT_m','cruce_PT'],kind='stable').reset_index(drop=True)
    seq=d.groupby('rid').cumcount()+1
    d.insert(0,'id_medicion',[f'R{int(r):02d}-T{int(n):03d}' for r,n in zip(d.rid,seq)])
    assert d.id_medicion.is_unique and d.H_m.eq(d.H_DEM_terraza_m).all()
    return refs,d,qc

def capas_revision(base,apoyo,catalogo):
    empty=lambda:gpd.GeoDataFrame({'rid':pd.Series(dtype=int)},geometry=gpd.GeoSeries([],crs=32718),crs=32718)
    capas={k:empty() for k in ['bordes','prolongaciones','transectos']}
    # Sólo cargar capas efectivamente dibujadas: nunca sustituir un borde ausente por una isohipsa.
    gpkg=Path(base)/'Revision_B_20261016/Revision_B_EDITAR_EN_QGIS.gpkg'
    if gpkg.exists():
        import pyogrio
        layers={x[0] for x in pyogrio.list_layers(gpkg)}
        for k in capas:
            if k in layers:
                g=gpd.read_file(gpkg,layer=k).to_crs(32718)
                if len(g):capas[k]=g
    jm=catalogo.copy();jm['x']=jm.este_utm18s_m;jm['y']=jm.norte_utm18s_m
    capas['jm']=jm;capas['iniciales']=gpd.read_file(Path(apoyo)/'Fuentes/AnalisisK.shp')
    return capas

def adaptar_figura(refs,res):
    r=res.copy();r['K_George_m02_anio']=r.K_m02_anio
    r['PT_elegible']=r.K_m02_anio.gt(0)&np.isfinite(r.K_m02_anio)
    r['principal_operativo']=False # ninguna elección automática de PT
    return refs.copy(),r

def tramo_ajuste(p,s0,s1):
    s=p.distancia_m.to_numpy();inner=s[(s>s0+1e-8)&(s<s1-1e-8)]
    ss=np.r_[s0,inner,s1]
    real=int(((s>=s0-1e-8)&(s<=s1+1e-8)).sum())
    x=np.interp(ss,s,p.chi_m);z=np.interp(ss,s,p.elevacion_m)
    return ss,x,z,real

def tasas_y_comparaciones(res,perfiles,catalogo,tau_KP,mar):
    cat=catalogo.set_index('terrace_id');rows=[]
    kp=tau_KP.set_index('id_medicion') if len(tau_KP) else pd.DataFrame()
    for _,v in res.iterrows():
        if not np.isfinite(v.K_m02_anio) or v.K_m02_anio<=0:continue
        p=perfiles[int(v.rid)];s,x,z,n=tramo_ajuste(p,v.s_B_m,v.s_PT_m)
        fit=rb.linregress(x,z) if n>=5 and np.ptp(x)>0 else None
        ksn=fit.slope*rb.A0**rb.THETA if fit is not None and fit.slope>0 else np.nan
        tau=(x-v.chi_B_m)/(v.K_m02_anio*rb.A0**rb.THETA)/1000
        u_tau=rb.linregress(tau,z).slope if fit is not None and fit.slope>0 else np.nan
        u_k=v.K_m02_anio*ksn*1000
        if np.isfinite(u_k):assert np.isclose(u_k,u_tau,rtol=1e-10,atol=1e-12)
        marker=cat.loc[v.terraza_id];sea=mar[v.mis]
        H=float(marker.elevacion_m);eH=float(marker.error_elevacion_m)
        T=float(v.T_ka);eT=float(v.error_T_ka)
        uterr=(H-sea['nivel_m'])/T
        eterr=np.sqrt((eH/T)**2+(sea['error_m']/T)**2+((H-sea['nivel_m'])*eT/T**2)**2)
        ujm=float(marker.tasa_alzamiento_m_ka);ejm=float(marker.error_tasa_alzamiento_m_ka)
        kjm=(ujm/1000)/ksn if ksn>0 and ujm>0 else np.nan
        kpend=(uterr/1000)/ksn if ksn>0 and uterr>0 else np.nan
        # Latitud del outlet de análisis original, sin cambiar al outlet costero recuperado.
        original=p.loc[~p.tramo_recuperado_aguas_abajo.astype(bool)] if 'tramo_recuperado_aguas_abajo' in p else p
        if original.empty:original=p
        pos=original.iloc[0];_,lat=LL.transform(pos.x_utm,pos.y_utm)
        trow=kp.loc[v.id_medicion] if len(kp) and v.id_medicion in kp.index else None
        dt=float(trow.chi_kp_m-v.chi_B_m) if trow is not None and trow.estado=='utilizable' else np.nan
        tk=float(trow.tau_KP_ka) if trow is not None else np.nan
        lim=ct.limites_grant()
        rows.append({**v.to_dict(),'latitud_rio':lat,'latitud_sur':-lat,'latitud_JM_sur':-float(marker.latitud_wgs84),
            'ksn_B_PT':ksn,'R2_B_PT':fit.rvalue**2 if fit else np.nan,'n_nodos_reales_ajuste':n,
            'U_aparente_Kksn_mm_anio':u_k,'U_aparente_tau_mm_anio':u_tau,
            'error_U_aparente_solo_K':u_k*v.error_K_propagado_aprox/v.K_m02_anio,
            'U_terraza_publicada_mm_anio':uterr,'error_U_terraza_mm_anio':eterr,
            'U_JM_publicado_mm_anio':ujm,'error_U_JM_publicado_mm_anio':ejm,
            'H_publicada_para_U_m':H,'nivel_mar_para_U_m':sea['nivel_m'],
            'U_DEM_provisional_mm_anio':(v.H_m-sea['nivel_m'])/T,
            'K_pendiente_terraza_m02_anio':kpend,'K_derivado_U_JM_m02_anio':kjm,
            'tau_KP_ka':tk,'error_tau_solo_K_ka':tk*v.error_K_propagado_aprox/v.K_m02_anio,
            'tau_KP_Grant_min_ka':dt/(lim['central_max']*rb.A0**rb.THETA)/1000,
            'tau_KP_Grant_max_ka':dt/(lim['central_min']*rb.A0**rb.THETA)/1000,
            'tau_KP_U_JM_ka':dt/(kjm*rb.A0**rb.THETA)/1000 if kjm>0 else np.nan,
            'U_Grant_min_mm_anio':lim['central_min']*ksn*1000,
            'U_Grant_max_mm_anio':lim['central_max']*ksn*1000,
            'nota_U':'dz/dtau = K*ksn; tasa aparente de incisión, igual a U sólo bajo equilibrio',
            'nota_K_JM':'K derivado aquí de U publicado J-M y ksn de este mismo tramo; no K publicado por J-M'})
    return pd.DataFrame(rows)

def figuras_perfil_K_tau(datos,perfiles,out):
    out=Path(out);out.mkdir(exist_ok=True);files=[]
    for rid,g in datos.groupby('rid'):
        p=perfiles[int(rid)];fig,axes=plt.subplots(2,1,figsize=(11,7),sharex=True,layout='constrained')
        for i,(bid,q) in enumerate(g.groupby('id_B_corto',sort=False)):
            color=ct.COLORS[i%len(ct.COLORS)]
            for _,v in q.iterrows():
                tail=p[p.chi_m>=v.chi_B_m];ss=(tail.distancia_m-v.s_B_m)/1000
                tt=(tail.chi_m-v.chi_B_m)/(v.K_m02_anio*rb.A0**rb.THETA)/1000
                axes[0].plot(ss,np.full(len(ss),v.K_m02_anio),color=color,alpha=.22,lw=.8)
                axes[1].plot(ss,tt,color=color,alpha=.25,lw=.8)
            axes[0].plot([],[],color=color,label=f'{bid} · MIS {q.mis.iloc[0]} · {len(q)} escenarios')
        axes[0].set_ylabel('K (m⁰·²/año)');axes[0].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        axes[1].set_ylabel('τ relativo a B (ka)');axes[1].set_xlabel('Distancia aguas arriba desde cada B (km)')
        for ax in axes:ax.grid(alpha=.2);ax.set_xlim(left=0);ax.set_ylim(bottom=0)
        axes[0].legend(fontsize=8);axes[0].set_title('Cada escenario supone K uniforme; no se infiere un K espacial nuevo')
        fig.suptitle(f'R{int(rid):02d} · {g.rio.iloc[0]} · K y τ a lo largo del cauce')
        f=out/f'R{int(rid):02d}_K_tau_distancia.png';fig.savefig(f,dpi=165,facecolor='white');plt.close(fig);files.append(f)
    return files

def comparacion_norte_sur(d,catalogo,fallas,out,mostrar):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    out=Path(out);out.mkdir(exist_ok=True)
    specs=[('K','K (m⁰·²/año)',[('K_m02_anio','K por B–PT','propio'),('K_derivado_U_JM_m02_anio','K derivado de U J-M','K_JM'),('K_pendiente_terraza_m02_anio','K por pendiente y U terraza','K_pendiente')]),
           ('Tau','τ en KP (ka)',[('tau_KP_ka','τ con K propios','propio'),('tau_KP_U_JM_ka','τ con K derivado de U J-M','K_JM')]),
           ('U','Tasa (mm/año)',[('U_aparente_tau_mm_anio','U aparente = dz/dτ = K×ksn','U_aparente'),('U_terraza_publicada_mm_anio','U de terraza (altura publicada)','U_terraza')])]
    interactive=make_subplots(rows=3,cols=1,shared_xaxes=True,vertical_spacing=.09,subplot_titles=[s[1] for s in specs])
    paths=[];summary=[]
    for panel,(name,label,series) in enumerate(specs,1):
        fig,ax=plt.subplots(figsize=(12,5),layout='constrained')
        for si,(column,title,option) in enumerate(series):
            if not mostrar.get(option,True):continue
            for mis,group in d.groupby('mis'):
                q=group[np.isfinite(group[column])].copy()
                if q.empty:continue
                color=('#7e42a7' if mis=='5e' else '#158663') if si==0 else (('#d17a12' if mis=='5e' else '#b84a16') if si==1 else '#656b74')
                marker=['o','^','s'][si];lab=title+f' · MIS {mis}'
                if mostrar.get('escenarios',True):
                    ax.scatter(q.latitud_sur,q[column],c=color,s=17,marker=marker,alpha=.4,zorder=3)
                stats=q.groupby('rid').agg(x=('latitud_sur','first'),mediana=(column,'median'),minimo=(column,'min'),maximo=(column,'max'),n=(column,'count')).sort_values('x')
                for rid,v in stats.iterrows():summary.append({'variable':column,'rid':rid,'MIS':mis,**v.to_dict()})
                x=stats.x.to_numpy();y=stats.mediana.to_numpy();low=stats.minimo.to_numpy();high=stats.maximo.to_numpy()
                ax.plot(x,y,color=color,marker=marker,ms=5,ls='-' if mostrar.get('lineas',True) else '',lw=1.2,label=lab)
                ax.fill_between(x,low,high,color=color,alpha=.10)
                interactive.add_trace(go.Scatter(x=q.latitud_sur.tolist(),y=q[column].tolist(),mode='markers',marker={'size':5,'color':color,'opacity':.4},name=lab+' · escenarios',legendgroup=lab,
                    customdata=q[['id_medicion','rio','id_B_corto','terraza_id']].to_numpy(),hovertemplate='%{customdata[1]} · %{customdata[0]}<br>B %{customdata[2]} · %{customdata[3]}<br>%{y:.5g}<extra></extra>'),row=panel,col=1)
                interactive.add_trace(go.Scatter(x=x.tolist(),y=y.tolist(),mode='lines+markers' if mostrar.get('lineas',True) else 'markers',line={'color':color},name=lab+' · mediana',legendgroup=lab),row=panel,col=1)
                if mostrar.get('errores_publicados',True) and si==0:
                    errcol={'K':'error_K_propagado_aprox','Tau':'error_tau_solo_K_ka','U':'error_U_aparente_solo_K'}[name]
                    qq=q[np.isfinite(q[errcol])];ax.errorbar(qq.latitud_sur,qq[column],yerr=qq[errcol],fmt='none',ecolor=color,alpha=.18,lw=.6)
        if mostrar.get('Grant',True):
            if name=='K':
                lim=ct.limites_grant();ax.axhspan(lim['outer_min'],lim['outer_max'],color='#658a7c',alpha=.10,label='Grant: envolvente con errores')
                ax.axhspan(lim['central_min'],lim['central_max'],color='#658a7c',alpha=.25,label='Grant: rango central')
                for low,high,opacity in [(lim['outer_min'],lim['outer_max'],.10),(lim['central_min'],lim['central_max'],.25)]:interactive.add_hrect(y0=low,y1=high,fillcolor='#658a7c',opacity=opacity,line_width=0,row=panel,col=1)
            else:
                lower,upper=('tau_KP_Grant_min_ka','tau_KP_Grant_max_ka') if name=='Tau' else ('U_Grant_min_mm_anio','U_Grant_max_mm_anio')
                q=d[np.isfinite(d[lower])&np.isfinite(d[upper])]
                for mis,qm in q.groupby('mis'):
                    g=qm.groupby('rid').agg(x=('latitud_sur','first'),low=(lower,'min'),high=(upper,'max')).sort_values('x')
                    ax.fill_between(g.x,g.low,g.high,color='#658a7c',alpha=.18,label='Referencia Grant · MIS '+mis)
                    interactive.add_trace(go.Scatter(x=g.x.tolist(),y=g.low.tolist(),mode='lines',line={'width':0},showlegend=False,legendgroup='Grant'+mis),row=panel,col=1)
                    interactive.add_trace(go.Scatter(x=g.x.tolist(),y=g.high.tolist(),mode='lines',line={'width':0},fill='tonexty',fillcolor='rgba(101,138,124,.18)',name='Referencia Grant · MIS '+mis,legendgroup='Grant'+mis),row=panel,col=1)
        if name=='U' and mostrar.get('JM_U_publicado',True):
            cat=catalogo[catalogo.terrace_id.isin(d.terraza_id)].drop_duplicates('terrace_id')
            for mis,q in cat.groupby('mis_asignado'):
                x=-q.latitud_wgs84.to_numpy();y=q.tasa_alzamiento_m_ka.to_numpy();e=q.error_tasa_alzamiento_m_ka.to_numpy()
                ax.errorbar(x,y,yerr=e if mostrar.get('errores_publicados',True) else None,fmt='*',color='#222',ms=8,label='U publicado J-M · '+mis,alpha=.75)
                interactive.add_trace(go.Scatter(x=x.tolist(),y=y.tolist(),mode='markers',marker={'symbol':'star','size':9,'color':'#222'},error_y={'type':'data','array':e.tolist(),'visible':mostrar.get('errores_publicados',True)},name='U publicado J-M · '+mis,customdata=q.terrace_id,hovertemplate='%{customdata}<br>U %{y:.3f} mm/año<extra></extra>'),row=panel,col=1)
        for system in ['Carranza','Pelluhue']:
            f=fallas[fallas.falla.astype(str).str.replace('*','',regex=False).str.strip().eq(system)]
            if f.empty:continue
            coords=np.concatenate([np.array(line.coords) for g in f.geometry for line in (list(g.geoms) if g.geom_type=='MultiLineString' else [g])]);west=coords[np.argmin(coords[:,0])];_,lat=LL.transform(*west)
            ax.axvline(-lat,ls='--',color='#555',lw=.8);ax.text(-lat,.98,system,rotation=90,va='top',ha='right',transform=ax.get_xaxis_transform(),fontsize=8)
            interactive.add_vline(x=-lat,line_dash='dash',line_color='#777',row=panel,col=1)
        ax.set_xlabel('Latitud (°S) · norte → sur');ax.set_ylabel(label);ax.grid(alpha=.2)
        if name=='K':ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        ax.set_title(name+' de norte a sur · línea = mediana descriptiva; franja = rango entre escenarios')
        ax.legend(fontsize=7,ncol=2,loc='best')
        path=out/f'{name}_norte_sur.png';fig.savefig(path,dpi=175,facecolor='white');plt.close(fig);paths.append(path)
    interactive.update_layout(template='plotly_white',height=1400,title='K, τ y tasas aparentes de norte a sur',legend={'orientation':'h','y':-.10},margin={'b':220})
    interactive.update_xaxes(title_text='Latitud (°S) · norte → sur',row=3,col=1)
    html=out/'Comparacion_K_tau_U_norte_sur.html';interactive.write_html(html,include_plotlyjs=True,config={'scrollZoom':True,'displaylogo':False})
    pd.DataFrame(summary).to_csv(out/'Resumen_por_rio_MIS.csv',index=False,encoding='utf-8-sig')
    return paths,html

def resumen_png(datos,out):
    out=Path(out);out.mkdir(exist_ok=True)
    rows=[]
    for (rid,bid),q in datos.groupby(['rid','id_B_corto'],sort=False):
        v=q.iloc[0];f=lambda col:f'{q[col].min():.3e} a {q[col].max():.3e}' if np.isfinite(q[col]).any() else 'Sin KP'
        rows.append([f'R{int(rid):02d}',bid,'5'+v.mis[1:],len(q),f('K_m02_anio'),f('tau_KP_ka'),f('U_aparente_tau_mm_anio'),f'{v.U_terraza_publicada_mm_anio:.3f}'])
    fig,ax=plt.subplots(figsize=(17,max(4,1+.32*len(rows))));ax.axis('off')
    table=ax.table(cellText=rows,colLabels=['Río','B','MIS','N','K (m⁰·²/año)','τ KP (ka)','U aparente (mm/año)','U terraza (mm/año)'],loc='center',cellLoc='center',colWidths=[.05,.04,.04,.035,.23,.20,.23,.13])
    table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.5)
    for (i,j),cell in table.get_celld().items():
        cell.set_edgecolor('#d2dee6');cell.set_facecolor('#214d63' if i==0 else '#eef3f6' if i%2==0 else 'white');cell.get_text().set_color('white' if i==0 else '#18394a')
    fig.suptitle('Resumen por B: todos los escenarios, sin elegir K definitivo',fontsize=13)
    path=out/'Resumen_K_tau_U_por_B.png';fig.savefig(path,dpi=175,bbox_inches='tight',facecolor='white');plt.close(fig);return path
