"""K con B sólo planimétrico y H del DEM en la terraza J-M asociada."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import json,html
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter,MaxNLocator
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs
import escenarios_k_manuales as km
import revision_b_base as rb

PALETA=['#1764a0','#e17c13','#7e42a7','#158663','#c64359','#82602b']

def referencias_dem(casos,terrazas,nombres,dem):
    # Se conserva la altura publicada sólo como metadato de contraste.
    refs=km.asociar_alturas(casos,terrazas,nombres)
    refs=refs[refs.fuente_H.eq('JM')].copy().reset_index(drop=True)
    refs['H_publicada_JM_m']=refs.H_JM_m
    refs['H_DEM_terraza_m']=[km.muestra_dem(dem,(x,y)) for x,y in zip(refs.x_terraza,refs.y_terraza)]
    if not np.isfinite(refs.H_DEM_terraza_m).all():raise ValueError('Terraza sin cobertura DEM: no sustituir por la altura de B ni por H publicada.')
    refs['H_m']=refs.H_DEM_terraza_m
    refs['fuente_H']='DEM_en_terraza_JM'
    refs['error_H_m']=0. # DEM fijo, sin incertidumbre asignada.
    refs['caso_altura']=refs.caso_B+'_H_DEM_TERRAZA'
    refs['nota_H']='DEM original muestreado bilinealmente en XY de la terraza J-M más cercana del mismo MIS. Z de B no fija H.'
    refs['diferencia_DEM_menos_publicada_m']=refs.H_DEM_terraza_m-refs.H_publicada_JM_m
    refs['compatibilidad_vertical_interna']='H y perfil proceden del mismo DEM; datum absoluto sin confirmar'
    refs=refs.rename(columns={'H_dem_punto_m':'Z_DEM_punto_B_solo_diagnostico_m'})
    return refs

def figuras_distribucion(res,refs,out,mostrar=None,bins='auto',mostrar_promedio=True):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);mostrar=mostrar or {}
    d=res.copy().sort_values(['rid','caso_B','s_PT_m','cruce_PT'],kind='stable').reset_index(drop=True)
    seq=d.groupby('rid').cumcount()+1
    d.insert(0,'id_medicion',[f'R{int(r):02d}-T{int(n):03d}' for r,n in zip(d.rid,seq)])
    mask=d.K_m02_anio.gt(0)&np.isfinite(d.K_m02_anio)
    for field,choices in [('mis',{'5e':'MIS5e','5c':'MIS5c'}),('origen_B',{'anterior':'B_anteriores','nuevo 3D':'B_nuevos'}),
                           ('direccion',{'ascendente':'ascendentes','descendente':'descendentes'})]:
        mask &= d[field].isin([value for value,key in choices.items() if mostrar.get(key,True)])
    d.to_csv(out/'Todos_los_escenarios_v8_con_ID.csv',index=False,encoding='utf-8-sig')
    used=d[mask].copy();paths=[];parts=[];summaries=[];hist_data=[]
    for rid,rr in refs.groupby('rid',sort=True):
        g=used[used.rid.eq(rid)].copy().reset_index(drop=True)
        if g.empty:
            parts.append(f'<h2>{html.escape(rr.rio.iloc[0])}</h2><p>Sin cruces utilizables con los filtros actuales.</p>');continue
        mean=float(g.K_m02_anio.mean());above=int((g.K_m02_anio>mean).sum());below=int((g.K_m02_anio<mean).sum())
        used.loc[used.rid.eq(rid),'promedio_K_rio']=mean
        summaries.append({'rid':rid,'rio':g.rio.iloc[0],'n':len(g),'K_promedio':mean,'K_min':g.K_m02_anio.min(),
                         'K_max':g.K_m02_anio.max(),'n_sobre_promedio':above,'n_bajo_promedio':below})
        edges=np.histogram_bin_edges(g.K_m02_anio.to_numpy(),bins=bins)
        # Evitar más intervalos que mediciones en muestras pequeñas.
        if len(edges)>min(21,len(g)+1):edges=np.linspace(g.K_m02_anio.min(),g.K_m02_anio.max(),min(21,len(g)+1))
        order=g.id_B_corto.drop_duplicates().tolist();colors={b:PALETA[i%len(PALETA)] for i,b in enumerate(order)}
        fig,(ax,ah)=plt.subplots(1,2,figsize=(max(14,min(21,.2*len(g)+6)),7),
                                  gridspec_kw={'width_ratios':[2.25,1]},layout='constrained')
        interactive=make_subplots(rows=1,cols=2,column_widths=[.69,.31],horizontal_spacing=.09,
                                  subplot_titles=['K por ID de medición','Histograma de K'])
        positions={value:i for i,value in enumerate(g.id_medicion)}
        legends=[];groups=[];group_colors=[]
        for bid,q in g.groupby('id_B_corto',sort=False):
            color=colors[bid];v=q.iloc[0];origin='anterior' if v.origen_B=='anterior' else 'nuevo'
            legends.append(Line2D([],[],marker='o',ls='',color=color,label=f'{bid} · B {origin} · MIS {v.mis}'))
            groups.append(q.K_m02_anio.to_numpy());group_colors.append(color)
            for direction,group in q.groupby('direccion'):
                filled=direction=='ascendente';symbol='o' if v.mis=='5e' else 'D'
                ax.scatter([positions[i] for i in group.id_medicion],group.K_m02_anio,s=49,marker=symbol,
                           facecolors=color if filled else 'none',edgecolors=color,linewidths=1.2,zorder=3)
                custom=group[['id_B_corto','origen_B','mis','direccion','terraza_id','H_m','H_publicada_JM_m',
                              'T_ka','chi_B_m','chi_PT_m','delta_chi_m','cruce_PT']].to_numpy().tolist()
                interactive.add_trace(go.Scatter(x=group.id_medicion.tolist(),y=group.K_m02_anio.tolist(),mode='markers',
                    marker={'color':color,'size':9,'symbol':('circle' if v.mis=='5e' else 'diamond')+('' if filled else '-open'),
                            'line':{'color':color,'width':1.2}},name=f'{bid} · MIS {v.mis} · '+('ascendente' if filled else 'descendente'),
                    customdata=custom,legendgroup=bid,
                    hovertemplate='%{x}<br>K: %{y:.5e} m^0.2/año<br>B: %{customdata[0]} · %{customdata[1]}<br>MIS: %{customdata[2]} · %{customdata[3]}'
                        '<br>Terraza: %{customdata[4]}<br>H DEM en terraza: %{customdata[5]:.3f} m<br>H publicada (referencia): %{customdata[6]:.2f} m'
                        '<br>T MIS: %{customdata[7]} ka<br>χ(B): %{customdata[8]:.2f} m<br>χ(P_T): %{customdata[9]:.2f} m'
                        '<br>Δχ: %{customdata[10]:.2f} m<br>Cruce P_T: %{customdata[11]}<extra></extra>'),row=1,col=1)
            counts,_=np.histogram(q.K_m02_anio,bins=edges)
            for j,count in enumerate(counts):hist_data.append({'rid':rid,'B':bid,'bin_inferior':edges[j],'bin_superior':edges[j+1],'n':int(count)})
            interactive.add_trace(go.Bar(x=((edges[:-1]+edges[1:])/2).tolist(),y=counts.tolist(),width=np.diff(edges).tolist(),
                marker_color=color,name=bid,legendgroup=bid,showlegend=False,
                customdata=np.c_[edges[:-1],edges[1:]].tolist(),
                hovertemplate='K: %{customdata[0]:.3e} a %{customdata[1]:.3e}<br>Frecuencia: %{y}<extra>'+bid+'</extra>'),row=1,col=2)
        ah.hist(groups,bins=edges,stacked=True,color=group_colors,alpha=.8,edgecolor='white',linewidth=.5)
        if mostrar_promedio:
            ax.axhline(mean,color='#b11c32',ls='--',lw=1.7);ah.axvline(mean,color='#b11c32',ls='--',lw=1.7)
            legends.append(Line2D([],[],color='#b11c32',ls='--',label=f'Promedio = {mean:.3e}'))
            interactive.add_hline(y=mean,line_color='#b11c32',line_dash='dash',row=1,col=1)
            interactive.add_vline(x=mean,line_color='#b11c32',line_dash='dash',row=1,col=2)
        ax.set_xticks(range(len(g)),g.id_medicion,rotation=90,fontsize=7 if len(g)>40 else 8)
        ax.set_xlim(-.7,len(g)-.3);ax.set_ylim(bottom=0);ax.set_xlabel('ID de medición · orden sin significado espacial')
        ax.set_ylabel(r'$K$ (m$^{0.2}$ año$^{-1}$)');ah.set_xlabel(r'$K$ (m$^{0.2}$ año$^{-1}$)');ah.set_ylabel('Número de escenarios')
        for axis in [ax.yaxis,ah.xaxis]:
            formatter=ScalarFormatter(useMathText=True);formatter.set_powerlimits((0,0));axis.set_major_formatter(formatter)
        ah.yaxis.set_major_locator(MaxNLocator(integer=True));ax.grid(axis='y',alpha=.2);ah.grid(axis='y',alpha=.2)
        legends += [Line2D([],[],marker='o',ls='',color='#333',label='Círculo: MIS 5e'),
                    Line2D([],[],marker='D',ls='',color='#333',label='Rombo: MIS 5c'),
                    Line2D([],[],marker='o',ls='',mfc='#333',color='#333',label='Relleno: P_T ascendente'),
                    Line2D([],[],marker='o',ls='',mfc='none',color='#333',label='Vacío: P_T descendente')]
        fig.legend(handles=legends,loc='outside lower center',ncol=3,fontsize=8)
        fig.suptitle(f'R{int(rid):02d} · {g.rio.iloc[0]} — distribución de K\nH del DEM en la terraza J-M asociada · {len(g)} escenarios · promedio {mean:.3e}',fontsize=13)
        path=out/f'R{int(rid):02d}_K_dispersion_histograma.png';fig.savefig(path,dpi=170,facecolor='white');plt.close(fig);paths.append(path)
        interactive.update_layout(template='plotly_white',height=720,barmode='stack',
            title={'text':f'{g.rio.iloc[0]} · {len(g)} escenarios · promedio {mean:.4e} m^0.2/año','font':{'size':17}},
            legend={'orientation':'h','y':-.35,'yanchor':'top','font':{'size':10},'itemclick':False,'itemdoubleclick':False},
            margin={'l':75,'r':25,'t':90,'b':230})
        interactive.update_xaxes(type='category',categoryorder='array',categoryarray=g.id_medicion.tolist(),tickangle=-90,
                                 title_text='ID de medición',row=1,col=1)
        interactive.update_yaxes(title_text='K (m^0.2/año)',rangemode='tozero',tickformat='.2e',row=1,col=1)
        interactive.update_xaxes(title_text='K (m^0.2/año)',tickformat='.2e',row=1,col=2)
        interactive.update_yaxes(title_text='Frecuencia',dtick=1,rangemode='tozero',row=1,col=2)
        block=interactive.to_html(include_plotlyjs=False,full_html=False,div_id=f'grafico-R{int(rid):02d}',
                                  config={'displaylogo':False,'responsive':True,'scrollZoom':True})
        details=' · '.join(f'{v.id_B_corto}: {v.terraza_id}, H_DEM={v.H_m:.2f} m' for _,v in rr.drop_duplicates('id_B_corto').iterrows())
        parts.append(f'<section id="R{int(rid):02d}"><h2>{html.escape(g.rio.iloc[0])}</h2><p>{html.escape(details)}</p>'+block+
                     f'<p><a href="{path.name}">Abrir figura PNG</a> · Sobre el promedio: {above}; bajo: {below}. '
                     'El histograma cuenta los mismos puntos. Promedio descriptivo de los escenarios.</p></section>')
    if used.empty:raise ValueError('No hay K con los filtros seleccionados.')
    used.to_csv(out/'K_graficados_con_ID.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(summaries).to_csv(out/'Resumen_promedios.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(hist_data).to_csv(out/'Conteos_histogramas.csv',index=False,encoding='utf-8-sig')
    nav=' '.join(f'<a href="#R{int(r.rid):02d}">{html.escape(r.rio)}</a>' for r in pd.DataFrame(summaries).itertuples())
    doc='''<!doctype html><html lang="es"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>K — altura DEM de la terraza, dispersión e histogramas</title><style>body{font:15px system-ui;color:#183247;background:white;max-width:1600px;margin:25px auto;padding:0 20px}section{margin:35px 0}nav{display:flex;gap:15px;flex-wrap:wrap}a{color:#175b8e}.nota{background:#eef4f8;padding:12px;border-radius:6px}</style>
<h1>K por río: dispersión e histograma</h1><p class="nota">B se utiliza sólo para fijar χ(B). H se obtiene del DEM en las coordenadas de la terraza J-M más cercana del mismo MIS, y define los cruces P_T. La altura del clic B no se usa para H.</p>
<p>Color = B específico (A: anterior; N: nuevo). Círculo = MIS 5e; rombo = MIS 5c. Relleno = cruce P_T ascendente; vacío = descendente. Línea roja = promedio aritmético. Señala un punto para ver sus datos; puedes ampliar la gráfica. Histogramas apilados por B, sin tablas.</p><nav>'''+nav+'</nav><script>'+get_plotlyjs()+'</script>'+''.join(parts)+'</html>'
    target=out/'ABRIR_K_TERRAZA_DEM.html';target.write_text(doc,encoding='utf8')
    return used,pd.DataFrame(summaries),paths,target
