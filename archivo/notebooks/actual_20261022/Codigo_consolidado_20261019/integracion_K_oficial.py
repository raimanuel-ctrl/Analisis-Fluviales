"""Aplica la revisión manual de B al cálculo 03, conservando roles y tres P_T."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROLES=('oficial','comparacion','prueba')
ETIQUETAS=(('primero','P1'),('intermedio_real','PM'),('ultimo','P3'))

def aplicar_decisiones(casos, revision):
    revision=Path(revision)
    d=pd.read_csv(revision/'Decisiones_B_20261019.csv',dtype={'mis':str},float_precision='round_trip').fillna('')
    if not d.id_B_corto.is_unique or not casos.id_B_corto.is_unique:
        raise ValueError('La revisión requiere un identificador único por posición B.')
    if set(d.id_B_corto)!=set(casos.id_B_corto):
        raise ValueError('Los B de entrada difieren de los revisados: revisar los nuevos puntos antes de asignarles rol.')
    if not set(d.rol_resultado).issubset(ROLES) or d.decision_B.eq('pendiente').any():
        raise ValueError('Hay roles inválidos o B pendientes en la revisión.')
    source=casos.set_index('id_B_corto',drop=False);rows=[]
    for _,decision in d.iterrows():
        b=source.loc[decision.id_B_corto].copy()
        for field in ['x_B','y_B','chi_B_m','s_B_m','z_B_m']:
            if not np.isclose(float(b[field]),float(decision[field]),rtol=0,atol=1e-6):
                raise ValueError(f'{decision.id_B_corto}: cambió {field}; repetir la revisión geométrica.')
        original=b.copy()
        b['mis_fuente']=str(original['mis']);b['caso_B_fuente']=original['caso_B']
        for field in ['mis','caso_B','rol_resultado','decision_B','estado_calibracion','comentario_usuario','fecha_decision']:
            b[field]=decision[field]
        b['B_oficial']=b['rol_resultado']=='oficial'
        b['revisado']=b['B_oficial']
        b['estado_B']={'oficial':'B oficial por selección manual del usuario',
          'comparacion':'alternativa conservada para comparación',
          'prueba':'posición B conservada; asignación MIS en prueba'}[b['rol_resultado']]
        b['clasificacion_MIS']='MIS según revisión manual explícita; se conserva MIS fuente por separado'
        rows.append(b.to_dict())
        if str(original['mis'])!=str(decision['mis']):
            original['mis_fuente']=str(original['mis']);original['caso_B_fuente']=original['caso_B']
            original['rol_resultado']='comparacion';original['decision_B']='comparacion'
            original['estado_calibracion']='historico';original['B_oficial']=False
            original['comentario_usuario']='Asignación MIS anterior conservada sólo como comparación histórica.'
            rows.append(original.to_dict())
    out=pd.DataFrame(rows)
    if not out.caso_B.is_unique:raise ValueError('Caso B/MIS duplicado tras aplicar decisiones.')
    return out,d

def seleccionar_tres(resultados):
    usable=resultados[np.isfinite(resultados.K_m02_anio)&resultados.K_m02_anio.gt(0)].copy()
    rows=[]
    for case,g in usable.groupby('caso_B',sort=False):
        g=g.sort_values(['s_PT_m','cruce_PT'],kind='stable').reset_index(drop=True)
        first,last=g.iloc[0],g.iloc[-1];target=(first.s_PT_m+last.s_PT_m)/2
        mid=g.iloc[int((g.s_PT_m-target).abs().argmin())]
        distinct=len({float(first.s_PT_m),float(mid.s_PT_m),float(last.s_PT_m)})
        for order,((label,suffix),v) in enumerate(zip(ETIQUETAS,[first,mid,last])):
            r=v.to_dict();T=float(v.T_ka);dt=float(v.error_T_ka)
            if not T>dt>=0:raise ValueError('Edad/rango cronológico inválido.')
            r.update(escenario_PT=label,orden_PT=order,
              id_medicion_cruce=v.id_medicion,escenario_K_cruce=v.escenario_K,
              id_medicion=f'R{int(v.rid):02d}-{v.id_B_corto}-{v.mis}-{suffix}',
              escenario_K=f'{v.escenario_K}_{label}',id_escenario_revision=f'R{int(v.rid):02d}_{v.id_B_corto}_{v.mis}_{label}',
              n_interceptos_reales=len(g),n_PT_distintos_retenidos=distinct,
              s_media_extremos_m=target,diferencia_respecto_distancia_media_m=v.s_PT_m-target,
              distancia_B_PT_m=v.s_PT_m-v.s_B_m,
              K_T_mas=v.K_m02_anio*T/(T+dt),K_T_menos=v.K_m02_anio*T/(T-dt),
              K_rango_geometrico_min=first.K_m02_anio,K_rango_geometrico_max=last.K_m02_anio,
              K_envolvente_geometria_edad_min=first.K_m02_anio*T/(T+dt),
              K_envolvente_geometria_edad_max=last.K_m02_anio*T/(T-dt),
              coincide_con_otro_escenario=distinct<3,
              estado='tres escenarios P_T; rol '+str(v.rol_resultado)+'; sin K único',
              nota='Tres etiquetas de P_T; las coincidencias no son observaciones independientes. Sin K promedio entre MIS.')
            rows.append(r)
    selected=pd.DataFrame(rows)
    if not selected.id_medicion.is_unique or not selected.escenario_K.is_unique:
        raise ValueError('Identificadores de escenarios no únicos.')
    return selected.sort_values(['rid','mis','caso_B','orden_PT'],kind='stable').reset_index(drop=True)

def verificar_revision(selected,revision,perfiles,folder):
    revision,folder=map(Path,[revision,folder]);folder.mkdir(exist_ok=True,parents=True)
    expected=pd.read_csv(revision/'Tres_PT_todos_los_rios_20261019.csv',dtype={'mis':str},float_precision='round_trip')
    expected=expected[expected.rid.isin(selected.rid.unique())]
    cols=['caso_B','escenario_PT','id_B_corto','mis','rol_resultado','terraza_id','x_B','y_B','s_B_m','chi_B_m',
          's_PT_m','chi_PT_m','H_m','T_ka','error_T_ka','K_m02_anio']
    joined=selected[cols].merge(expected[cols],on=['caso_B','escenario_PT'],suffixes=('_calculado','_revision'),validate='one_to_one',how='outer',indicator=True)
    if not joined['_merge'].eq('both').all():raise ValueError('Faltan o sobran casos respecto de la revisión guardada.')
    for col in ['id_B_corto','mis','rol_resultado','terraza_id']:
        if not joined[col+'_calculado'].eq(joined[col+'_revision']).all():raise ValueError('Cambió '+col)
    for col in ['x_B','y_B','s_B_m','chi_B_m','s_PT_m','chi_PT_m','H_m','T_ka','error_T_ka','K_m02_anio']:
        a=joined[col+'_calculado'];b=joined[col+'_revision']
        tol=0 if col=='K_m02_anio' else 1e-7
        if not np.allclose(a,b,rtol=1e-12,atol=tol):raise ValueError('El recálculo no reproduce '+col)
        joined['diferencia_'+col]=a-b
    residual=[]
    for v in selected.itertuples():
        p=perfiles[int(v.rid)];residual.append(abs(np.interp(v.s_PT_m,p.distancia_m,p.elevacion_m)-v.H_m))
    if max(residual)>=1e-6:raise ValueError('Algún P_T no coincide con H.')
    if not selected.groupby('caso_B').size().eq(3).all():raise ValueError('Se requieren tres etiquetas por caso.')
    joined.to_csv(folder/'Control_reproduccion_revision.csv',index=False,encoding='utf-8-sig')
    report={'n_escenarios':len(selected),'n_oficiales':int(selected.rol_resultado.eq('oficial').sum()),
      'n_B_MIS_oficiales':int(selected[selected.rol_resultado.eq('oficial')].caso_B.nunique()),
      'n_comparacion':int(selected.rol_resultado.eq('comparacion').sum()),'n_prueba':int(selected.rol_resultado.eq('prueba').sum()),
      'max_error_H_PT_m':float(max(residual)),
      'max_error_relativo_K':float(np.max(np.abs(joined.K_m02_anio_calculado/joined.K_m02_anio_revision-1))),
      'IDs_unicos':True,'promedios_entre_MIS':False,'K_unico_adoptado':False}
    (folder/'Control_seleccion_oficial.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    return report

def exportar_roles(casos,refs,raw,selected,revision,out):
    revision,out=map(Path,[revision,out]);archive=out/'Archivo_escenarios';archive.mkdir(exist_ok=True)
    casos.to_csv(out/'B_todos_roles.csv',index=False,encoding='utf-8-sig')
    raw.to_csv(archive/'K_todos_los_interceptos_todos_roles.csv',index=False,encoding='utf-8-sig')
    refs.to_csv(archive/'Referencias_todos_roles.csv',index=False,encoding='utf-8-sig')
    selected.to_csv(archive/'K_tres_PT_todos_roles.csv',index=False,encoding='utf-8-sig')
    for role in ROLES:
        root=out if role=='oficial' else out/('Comparacion' if role=='comparacion' else 'Prueba')
        root.mkdir(exist_ok=True)
        subset=selected[selected.rol_resultado.eq(role)]
        subset.to_csv(root/'K_tres_PT.csv',index=False,encoding='utf-8-sig')
        raw[raw.rol_resultado.eq(role)].to_csv(root/'K_todos_interceptos.csv',index=False,encoding='utf-8-sig')
        casos[casos.rol_resultado.eq(role)].to_csv(root/'B_utilizados.csv',index=False,encoding='utf-8-sig')
        refs[refs.rol_resultado.eq(role)].to_csv(root/'Referencias_H_DEM_terraza.csv',index=False,encoding='utf-8-sig')
    official=selected[selected.rol_resultado.eq('oficial')]
    # Compatibilidad: este nombre contiene los tres escenarios oficiales usados por Tau/.
    official.to_csv(out/'K_todos_los_cruces.csv',index=False,encoding='utf-8-sig')
    official.to_csv(out/'K_oficiales_tres_PT.csv',index=False,encoding='utf-8-sig')
    snapshot=out/'Fuentes_revision';snapshot.mkdir(exist_ok=True)
    for name in ['Decisiones_B_20261019.csv','Tres_PT_todos_los_rios_20261019.csv','K_oficiales_tres_PT_20261019.csv']:
        import shutil
        shutil.copy2(revision/name,snapshot/name)
    (out/'LEEME_alcance_archivos.txt').write_text(
      'Resultados principales: sólo B oficiales, tres etiquetas P_T por B/MIS.\n'
      'K_todos_los_cruces.csv conserva el nombre compatible y contiene los mismos escenarios que Tau/.\n'
      'Todos los interceptos oficiales están en K_todos_interceptos.csv.\n'
      'Archivo_escenarios/ conserva todos los interceptos y los tres P_T de todos los roles.\n'
      'Comparacion/ y Prueba/ tienen sus propios B, K y productos; no alimentan los resultados oficiales.\n'
      'Las coincidencias geométricas entre etiquetas no son mediciones independientes.\n',encoding='utf-8')

def figuras_K(datos,out,mostrar,grant):
    """Dispersión e histograma por río/MIS; ningún promedio entre MIS ni K único."""
    import matplotlib.pyplot as plt
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    from plotly.offline import get_plotlyjs
    from html import escape
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=datos.copy();mask=d.K_m02_anio.gt(0)
    for field,mapping in [('mis',{'5e':'MIS5e','5c':'MIS5c'}),('origen_B',{'anterior':'B_anteriores','nuevo 3D':'B_nuevos'}),('direccion',{'ascendente':'ascendentes','descendente':'descendentes'})]:
        mask &= d[field].isin([k for k,v in mapping.items() if mostrar.get(v,True)])
    used=d[mask];figures=[];blocks=[];counts=[]
    palette=['#1764a0','#b87517','#8250a0','#21845b']
    for (rid,mis),g in used.groupby(['rid','mis'],sort=True):
        fig,(ax,ah)=plt.subplots(1,2,figsize=(13,5.8),layout='constrained')
        inter=make_subplots(rows=1,cols=2,subplot_titles=['Tres etiquetas de P_T por B','Histograma de etiquetas de escenario'])
        ordered=g.sort_values(['caso_B','orden_PT']);positions={v:i for i,v in enumerate(ordered.id_medicion)}
        values=g.K_m02_anio.to_numpy()
        if np.ptp(values)==0:
            center=float(values[0]);half=max(abs(center)*.05,1e-15)
            bins=np.array([center-half,center+half])
        else:
            bins=np.histogram_bin_edges(values,bins=min(8,max(1,len(g))))
        for j,(case,q) in enumerate(ordered.groupby('caso_B',sort=False)):
            color=palette[j%len(palette)];v=q.iloc[0]
            ax.scatter([positions[x] for x in q.id_medicion],q.K_m02_anio,color=color,s=45,label=v.id_B_corto)
            n,edges=np.histogram(q.K_m02_anio,bins=bins)
            ah.hist(q.K_m02_anio,bins=bins,alpha=.5,color=color,label=v.id_B_corto)
            for a,b,c in zip(edges[:-1],edges[1:],n):counts.append({'rid':rid,'mis':mis,'caso_B':case,'min_K':a,'max_K':b,'n_etiquetas':int(c)})
            inter.add_trace(go.Scatter(x=q.id_medicion.tolist(),y=q.K_m02_anio.tolist(),mode='markers',name=f'{v.id_B_corto} · {mis}',
                marker={'color':color,'size':9},customdata=q[['escenario_PT','n_PT_distintos_retenidos']].to_numpy(),
                hovertemplate='%{x}<br>K %{y:.6e}<br>%{customdata[0]}<br>Posiciones distintas del caso: %{customdata[1]}<extra></extra>'),row=1,col=1)
            inter.add_trace(go.Bar(x=((edges[1:]+edges[:-1])/2).tolist(),y=n.tolist(),width=np.diff(edges).tolist(),name=v.id_B_corto,marker_color=color,showlegend=False),row=1,col=2)
        for lo,hi,alpha in ([(grant['outer_min'],grant['outer_max'],.08),(grant['central_min'],grant['central_max'],.18)] if grant and np.isfinite(list(grant.values())).all() else []):
            ax.axhspan(lo,hi,color='#708a86',alpha=alpha);ah.axvspan(lo,hi,color='#708a86',alpha=alpha)
            inter.add_hrect(y0=lo,y1=hi,fillcolor='#708a86',opacity=alpha,line_width=0,row=1,col=1)
            inter.add_vrect(x0=lo,x1=hi,fillcolor='#708a86',opacity=alpha,line_width=0,row=1,col=2)
        ax.set_xticks(range(len(ordered)),ordered.id_medicion,rotation=70,ha='right',fontsize=8)
        ax.set_ylabel('K (m$^{0,2}$/año)');ax.set_xlabel('Etiqueta de escenario P_T');ah.set_xlabel('K (m$^{0,2}$/año)');ah.set_ylabel('Número de etiquetas')
        ax.legend(loc='lower left');ax.grid(axis='y',alpha=.2);ah.grid(axis='y',alpha=.2)
        title=f'{g.rio.iloc[0]} · MIS {mis} · {len(g)} etiquetas de P_T; sin promedio de K'
        fig.suptitle(title);target=out/f'R{int(rid):02d}_{mis}_K_dispersion_histograma.png'
        fig.savefig(target,dpi=165,facecolor='white');plt.close(fig);figures.append(target)
        inter.update_layout(template='plotly_white',height=640,title=title,barmode='stack',margin={'b':170})
        inter.update_xaxes(title_text='Escenario',tickangle=-70,row=1,col=1);inter.update_yaxes(title_text='K (m^0.2/año)',tickformat='.2e',row=1,col=1)
        inter.update_xaxes(title_text='K (m^0.2/año)',tickformat='.2e',row=1,col=2);inter.update_yaxes(title_text='Número de etiquetas',row=1,col=2)
        blocks.append(inter.to_html(full_html=False,include_plotlyjs=False,config={'displaylogo':False,'responsive':True}))
    used.to_csv(out/'K_graficados_con_ID.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(counts).to_csv(out/'Conteos_histogramas.csv',index=False,encoding='utf-8-sig')
    html=out/'ABRIR_K_TERRAZA_DEM.html'
    html.write_text('<!doctype html><meta charset="utf-8"><title>K por río y MIS</title><style>body{font:15px system-ui;margin:25px;color:#183247}section{margin:25px 0}</style><h1>K por río y MIS</h1><p>Sin promedios entre MIS. Los histogramas cuentan etiquetas de escenario; las coincidencias no son mediciones independientes. Franja: referencia Grant.</p><script>'+get_plotlyjs()+'</script>'+''.join(blocks),encoding='utf-8')
    return used,figures,html
