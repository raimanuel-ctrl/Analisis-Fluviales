"""Bandas automáticas N-S de ancho E-O fijo, relieve 3D y pendiente topográfica."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import json
import re
import warnings
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds
from scipy.ndimage import map_coordinates
import plotly.graph_objects as go
from plotly.subplots import make_subplots

COLORES_JM = {'5e':'#c10caf', '5c':'#138c45'}


def catalogo_terrazas(puntos):
    """Normaliza el catálogo indicado sin alterar alturas, MIS ni asociaciones."""
    t=puntos.copy()
    required={'terrace_id','este_utm18s_m','norte_utm18s_m','mis_asignado','elevacion_m',
              'error_elevacion_m','edad_referencia_ka','incertidumbre_edad_ka'}
    if not required.issubset(t):
        raise ValueError(f'Faltan campos en catálogo J-M: {required-set(t.columns)}')
    if t.terrace_id.isna().any() or t.terrace_id.duplicated().any():
        raise ValueError('El catálogo debe tener terrace_id únicos.')
    t['mis_vis']=t.mis_asignado.astype(str).str.lower().str.extract(r'5\s*([ec])',expand=False).map({'e':'5e','c':'5c'})
    if t.mis_vis.isna().any():
        raise ValueError('El catálogo contiene un MIS distinto de 5e/5c o sin identificar.')
    cols=['este_utm18s_m','norte_utm18s_m','elevacion_m','error_elevacion_m','edad_referencia_ka','incertidumbre_edad_ka']
    if not np.isfinite(t[cols].to_numpy(float)).all():
        raise ValueError('Catálogo con coordenadas, alturas o edades no finitas.')
    if (t[['error_elevacion_m','incertidumbre_edad_ka']]<0).any().any() or (t.edad_referencia_ka<=0).any():
        raise ValueError('Edades o incertidumbres inválidas en catálogo.')
    if 'crs_utm' in t and not t.crs_utm.astype(str).str.contains('32718').all():
        raise ValueError('Los puntos deben estar referidos a UTM18S EPSG:32718.')
    return t


def puntos_en_banda(banda, puntos, elevacion_visual_m=3.):
    t=catalogo_terrazas(puntos)
    t['este_desde_B_m']=t.este_utm18s_m-banda['x_B']
    t['norte_desde_B_m']=t.norte_utm18s_m-banda['y_B']
    t=t[t.este_desde_B_m.abs().le(banda['ancho_ew_m']/2) &
        t.norte_desde_B_m.abs().le(banda['largo_ns_m']/2)].copy()
    if t.empty:
        t['z_dem_m']=pd.Series(dtype=float)
    else:
        ex,ny=banda['east_m'],banda['north_m']
        cols=(t.este_desde_B_m.to_numpy()-ex[0])/(ex[1]-ex[0])
        rows=(t.norte_desde_B_m.to_numpy()-ny[0])/(ny[1]-ny[0])
        t['z_dem_m']=map_coordinates(banda['z_m'],[rows,cols],order=1,mode='constant',cval=np.nan,prefilter=False)
    t['z_simbolo_m']=t.z_dem_m+elevacion_visual_m
    t['simbolo_despegue_m']=elevacion_visual_m
    return t


def texto_terraza(row, incluir_dem=True):
    s=(f'{row.terrace_id} · MIS {row.mis_vis}<br>Altura publicada: {row.elevacion_m:.2f} ± '
       f'{row.error_elevacion_m:.2f} m<br>Edad MIS: {row.edad_referencia_ka:g} ± '
       f'{row.incertidumbre_edad_ka:g} ka<br>Este UTM: {row.este_utm18s_m:.1f} m'
       f'<br>Norte UTM: {row.norte_utm18s_m:.1f} m')
    if incluir_dem:
        s+=f'<br>Terreno DEM en XY: {row.z_dem_m:.2f} m<br>Marcador sobre DEM (+{row.simbolo_despegue_m:g} m visuales)'
    return s


def estadisticas_banda(z, coordenada, axis):
    """Estadísticas sin relleno de huecos. percentiles espaciales, no incertidumbre."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        p = np.nanpercentile(z, [0, 50, 90, 100], axis=axis)
    counts = np.isfinite(z).sum(axis=axis)
    return pd.DataFrame({'desplazamiento_m': coordenada, 'z_min_m': p[0],
        'z_mediana_m': p[1], 'z_p90_m': p[2], 'z_max_m': p[3],
        'n_validos': counts, 'fraccion_cobertura': counts/z.shape[axis]})


def preparar_banda(dem, x_B, y_B, semiancho_ew_m=500., semilargo_ns_m=1000.,
                   paso_m=None, max_nodos=350000):
    """Rectángulo N-S por B. Ancho = 2*semiancho E/O; longitud = 2*semilargo N/S.

    La pendiente se calcula con vecinos del DEM nativo antes de interpolar la malla 3D.
    No se usan alturas de terrazas para encontrar un borde ni se ajusta el datum.
    """
    if not np.isfinite([x_B, y_B, semiancho_ew_m, semilargo_ns_m]).all():
        raise ValueError('Centro y dimensiones deben ser finitos.')
    if semiancho_ew_m <= 0 or semilargo_ns_m <= 0:
        raise ValueError('Las dos semidimensiones deben ser positivas.')
    with rasterio.open(dem) as src:
        if src.crs != rasterio.crs.CRS.from_epsg(32718):
            raise ValueError('El DEM debe estar en EPSG:32718 para usar metros E/O y N/S.')
        if src.transform.b != 0 or src.transform.d != 0:
            raise ValueError('DEM rotado: reproyectar explícitamente antes de usar esta malla.')
        if paso_m is None:
            paso_m = max(abs(src.res[0]), abs(src.res[1]))
        if not np.isfinite(paso_m) or paso_m <= 0:
            raise ValueError('Paso de muestreo inválido.')
        nx = max(3, int(np.ceil(2*semiancho_ew_m/paso_m))+1)
        ny = max(3, int(np.ceil(2*semilargo_ns_m/paso_m))+1)
        # Número impar: la malla contiene el punto central exactamente.
        nx += nx % 2 == 0
        ny += ny % 2 == 0
        if nx*ny > max_nodos:
            raise ValueError(f'Malla de {nx*ny:,} nodos; aumentar paso_m o disminuir dimensiones.')
        east = np.linspace(-semiancho_ew_m, semiancho_ew_m, nx)
        north = np.linspace(-semilargo_ns_m, semilargo_ns_m, ny)
        xx, yy = np.meshgrid(x_B+east, y_B+north)
        halo = 3*max(src.res)
        win = from_bounds(x_B-semiancho_ew_m-halo, y_B-semilargo_ns_m-halo,
                          x_B+semiancho_ew_m+halo, y_B+semilargo_ns_m+halo,
                          src.transform).round_offsets().round_lengths()
        native = src.read(1, window=win, masked=True, boundless=True).astype(float).filled(np.nan)
        tr = src.window_transform(win)
        # Gradiente sobre el raster nativo, con halo que evita pendientes falsas en el límite.
        dzdy, dzdx = np.gradient(native, abs(src.res[1]), abs(src.res[0]))
        slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy)))
        inv = ~tr
        def sample(arr, x, y):
            cols = inv.a*x+inv.b*y+inv.c-.5
            rows = inv.d*x+inv.e*y+inv.f-.5
            return map_coordinates(arr, [rows, cols], order=1, mode='constant', cval=np.nan, prefilter=False)
        z = sample(native, xx, yy)
        slope_display = sample(slope, xx, yy)
        z_B = float(sample(native, np.array([x_B]), np.array([y_B]))[0])
        if not np.isfinite(z).any():
            raise ValueError('El rectángulo no contiene elevaciones válidas en el DEM.')
        return {'east_m': east, 'north_m': north, 'z_m': z, 'pendiente_grados': slope_display,
            'x_B': float(x_B), 'y_B': float(y_B), 'z_B_DEM_m': z_B,
            'ancho_ew_m': float(2*semiancho_ew_m), 'largo_ns_m': float(2*semilargo_ns_m),
            'paso_solicitado_m': float(paso_m), 'resolucion_dem_m': list(src.res),
            'cobertura': float(np.isfinite(z).mean()),
            'perfil_NS': estadisticas_banda(z, north, axis=1),
            'perfil_EW': estadisticas_banda(z, east, axis=0),
            'dem': str(Path(dem).resolve()),
            'nota': 'Pendiente topográfica del DEM; no ksn. Percentiles de banda espaciales, no IC.'}


def graficar_banda(banda, nombre, caso, perfil_rio=None, puntos_JM=None,
                   pendiente_max_color=45., exageracion_vertical=1., mostrar_H=None,
                   etiquetar_terrazas=True):
    if pendiente_max_color <= 0 or exageracion_vertical <= 0:
        raise ValueError('Máximo de pendiente y exageración deben ser positivos.')
    b = banda
    fig = make_subplots(rows=2, cols=2,
        specs=[[{'type':'scene', 'rowspan':2}, {'type':'xy'}], [None, {'type':'xy'}]],
        column_widths=[.64,.36], row_heights=[.5,.5], horizontal_spacing=.12, vertical_spacing=.17,
        subplot_titles=['Relieve 3D · color = pendiente topográfica',
                        f'Banda N–S · ancho O–E: {b["ancho_ew_m"]:g} m',
                        f'Banda O–E · largo N–S: {b["largo_ns_m"]:g} m'])
    east, north, z = b['east_m'], b['north_m'], b['z_m']
    xx, yy = np.meshgrid(east, north)
    # Una única escala común por defecto: 0° azul, 45° rojo; valores mayores saturan rojo.
    scale = [[0., '#1155c9'], [.35, '#39b9ef'], [.6, '#f1e8af'], [.8, '#f18c59'], [1., '#cb1828']]
    fig.add_trace(go.Surface(x=xx, y=yy, z=z, surfacecolor=b['pendiente_grados'],
        colorscale=scale, cmin=0, cmax=pendiente_max_color,
        colorbar={'title':{'text':'Pendiente<br>(°)','font':{'size':11}}, 'x':.58, 'y':.51, 'len':.68, 'thickness':15,
                  'tickvals':[0,pendiente_max_color/3,2*pendiente_max_color/3,pendiente_max_color],
                  'ticktext':['0',f'{pendiente_max_color/3:g}',f'{2*pendiente_max_color/3:g}',f'≥{pendiente_max_color:g}']},
        name='DEM', meta={'capa':'DEM'}, hovertemplate='Oeste/Este desde B: %{x:.1f} m<br>Sur/Norte desde B: %{y:.1f} m'
                                 '<br>Elevación DEM: %{z:.2f} m<extra>DEM</extra>',
        lighting={'ambient':.85,'diffuse':.65,'specular':.05,'roughness':.9},
        contours={'x':{'show':False},'y':{'show':False},'z':{'show':False}}),row=1,col=1)
    fig.add_trace(go.Scatter3d(x=[0],y=[0],z=[b['z_B_DEM_m']],mode='markers+text',
        text=['B inicial'],textposition='top center',marker={'size':7,'color':'black','symbol':'diamond'},
        name='B inicial por cota',meta={'capa':'B_inicial'},hovertemplate='B inicial<br>Z DEM: %{z:.2f} m<extra></extra>'),row=1,col=1)
    if perfil_rio is not None:
        ex = perfil_rio.x_utm.to_numpy()-b['x_B']
        ny = perfil_rio.y_utm.to_numpy()-b['y_B']
        valid = (abs(ex)<=b['ancho_ew_m']/2)&(abs(ny)<=b['largo_ns_m']/2)
        # Mantener NaN en salidas/reentradas: no conectar secciones externas al recorte.
        fig.add_trace(go.Scatter3d(x=np.where(valid,ex,np.nan),y=np.where(valid,ny,np.nan),
            z=np.where(valid,perfil_rio.elevacion_m,np.nan),mode='lines',connectgaps=False,
            line={'color':'#151515','width':6},name='Río actual',meta={'capa':'rio'},
            hovertemplate='Cauce actual<br>ΔE %{x:.0f} m · ΔN %{y:.0f} m<br>Z %{z:.2f} m<extra></extra>'),row=1,col=1)
    if puntos_JM is not None:
        pts=puntos_en_banda(b,puntos_JM)
        for mis,color in COLORES_JM.items():
            q=pts[pts.mis_vis.eq(mis) & pts.z_dem_m.notna()]
            fig.add_trace(go.Scatter3d(x=q.este_desde_B_m,y=q.norte_desde_B_m,z=q.z_simbolo_m,
                mode='markers+text' if etiquetar_terrazas else 'markers',text=q.terrace_id,
                textposition='top center',textfont={'size':10,'color':color},
                marker={'size':7,'color':color,'symbol':'diamond' if mis=='5e' else 'circle',
                        'line':{'color':'white','width':2}},
                name=f'Terrazas J-M MIS {mis}',legendgroup=f'JM_{mis}',
                meta={'capa':'terrazas_JM','mis':mis,'n_puntos':len(q)},
                hovertext=[texto_terraza(v) for _,v in q.iterrows()],hovertemplate='%{hovertext}<extra>J-M 2015</extra>'),row=1,col=1)
    # Eje central N-S y sección O-E de referencia, ambos sobre la superficie.
    midx, midy = len(east)//2, len(north)//2
    for ex, ny, zz, name in [(np.zeros(len(north)),north,z[:,midx],'Eje N–S por B'),
                              (east,np.zeros(len(east)),z[midy,:],'Sección O–E por B')]:
        fig.add_trace(go.Scatter3d(x=ex,y=ny,z=zz,mode='lines',line={'width':3,'color':'#fafafa'},
                                  name=name,meta={'capa':'eje_NS' if len(ex)==len(north) and np.all(ex==0) else 'eje_EW'},
                                  showlegend=False,hoverinfo='skip'),row=1,col=1)
    if mostrar_H is not None:
        fig.add_trace(go.Surface(x=xx[::8,::8],y=yy[::8,::8],z=np.full_like(z[::8,::8],mostrar_H),
            colorscale=[[0,'#a044af'],[1,'#a044af']],showscale=False,opacity=.17,name='Plano H provisional',
            meta={'capa':'plano_H','H_m':mostrar_H},hovertemplate=f'Plano H={mostrar_H:.2f} m · referencia vertical por comprobar<extra></extra>'),row=1,col=1)
    for row,key in [(1,'perfil_NS'),(2,'perfil_EW')]:
        q=b[key]
        xc=q.desplazamiento_m
        for col, style, name in [('z_min_m',{'color':'#9da8b3','width':1},'Mínimo'),
                                 ('z_max_m',{'color':'#9da8b3','width':1},'Máximo'),
                                 ('z_mediana_m',{'color':'#1b507c','width':2.5},'Mediana'),
                                 ('z_p90_m',{'color':'#a85c16','width':1.8,'dash':'dash'},'Percentil 90')]:
            fig.add_trace(go.Scatter(x=xc,y=q[col],mode='lines',line=style,name=name,
                legendgroup=name,showlegend=row==1,connectgaps=False,
                meta={'capa':'perfil_swath','direccion':'NS' if row==1 else 'EW','variable':col},
                customdata=q.fraccion_cobertura,
                hovertemplate='%{x:.0f} m · Z %{y:.2f} m<br>Cobertura %{customdata:.1%}<extra>'+name+'</extra>'),row=row,col=2)
        central_z=z[:,midx] if row==1 else z[midy,:]
        fig.add_trace(go.Scatter(x=xc,y=central_z,mode='lines',line={'color':'black','width':1.1,'dash':'dot'},
            name='Perfil central por B',legendgroup='central',showlegend=row==1,
            meta={'capa':'perfil_central','direccion':'NS' if row==1 else 'EW'}),row=row,col=2)
        fig.add_trace(go.Scatter(x=[0],y=[b['z_B_DEM_m']],mode='markers',name='B sobre el DEM',
            marker={'color':'black','size':7,'symbol':'diamond'},showlegend=False,
            meta={'capa':'B_perfil','direccion':'NS' if row==1 else 'EW'}),row=row,col=2)
        # add_vline inspecciona mal las trazas Surface en subplots mixtos en Plotly 7.1.
        fig.add_shape(type='line',x0=0,x1=0,y0=0,y1=1,
                      xref='x' if row==1 else 'x2',yref='y domain' if row==1 else 'y2 domain',
                      line={'width':1,'dash':'dot','color':'#5b6975'})
        fig.update_xaxes(title_text='Distancia N–S desde B (m)<br>Sur (−) · Norte (+)' if row==1 else
                         'Distancia O–E desde B (m)<br>Oeste (−) · Este (+)',row=row,col=2)
        fig.update_yaxes(title_text='Z DEM (m)',row=row,col=2)
    zspan=max(1.,float(np.nanmax(z)-np.nanmin(z)))
    scene={'xaxis':{'title':'Oeste (−) / Este (+), m'},
           'yaxis':{'title':'Sur (−) / Norte (+), m'},
           'zaxis':{'title':'Z DEM (m)'},
           'camera':{'eye':{'x':-.9,'y':-1.1,'z':.75}},'bgcolor':'#f5f7f9'}
    # Proporciones métricas explícitas: EV=1 preserva z/x e y/x; normalizar el lado
    # máximo evita que la normalización volumétrica de Plotly recorte bandas planas.
    scale=max(b['ancho_ew_m'],b['largo_ns_m'],zspan*exageracion_vertical)
    scene.update(aspectmode='manual',aspectratio={'x':b['ancho_ew_m']/scale,
        'y':b['largo_ns_m']/scale,'z':zspan*exageracion_vertical/scale})
    fig.update_layout(scene=scene,template='plotly_white',height=850,autosize=True,
        title={'text':f'{nombre} · {caso}<br><sup>Bandas de relieve: {b["ancho_ew_m"]:.0f} m O–E × '
                      f'{b["largo_ns_m"]:.0f} m S–N · azul = menor pendiente; rojo = mayor · '
                      f'EV={exageracion_vertical:g} · alturas DEM</sup>','x':.02},
        margin={'l':35,'r':35,'t':110,'b':105},
        legend={'orientation':'h','y':-.09,'x':.02,'font':{'size':11}},
        annotations=list(fig.layout.annotations)+[{'text':'Ctrl + clic: marcar punto · Arrastrar: rotar · Rueda: zoom · '
            'El color representa pendiente del terreno, no ksn ni una paleocosta detectada.',
            'xref':'paper','yref':'paper','x':0,'y':-.16,'showarrow':False,'font':{'size':11,'color':'#4b5e6d'}}])
    return fig


def mapa_contexto_terrazas(banda, puntos, perfil_rio=None, semilado_m=3000., etiquetar=True):
    """Mapa local: conserva el swath original y muestra puntos cercanos fuera de su recorte."""
    from matplotlib.colors import LightSource
    from PIL import Image
    import base64,io
    if semilado_m < max(banda['ancho_ew_m'],banda['largo_ns_m'])/2:
        raise ValueError('El mapa de contexto debe incluir toda la banda.')
    t=catalogo_terrazas(puntos)
    t['este_desde_B_m']=t.este_utm18s_m-banda['x_B']
    t['norte_desde_B_m']=t.norte_utm18s_m-banda['y_B']
    t=t[t.este_desde_B_m.abs().le(semilado_m) & t.norte_desde_B_m.abs().le(semilado_m)].copy()
    t['dentro_swath_3D']=t.este_desde_B_m.abs().le(banda['ancho_ew_m']/2) & t.norte_desde_B_m.abs().le(banda['largo_ns_m']/2)
    context=preparar_banda(banda['dem'],banda['x_B'],banda['y_B'],semilado_m,semilado_m,
                           paso_m=max(24.,banda['paso_solicitado_m']))
    z=context['z_m'][::-1]
    finite=np.isfinite(z)
    fill=np.where(finite,z,np.nanmedian(z))
    shade=LightSource(azdeg=315,altdeg=45).hillshade(fill,
        dx=context['east_m'][1]-context['east_m'][0],dy=context['north_m'][1]-context['north_m'][0])
    rgba=np.zeros(z.shape+(4,),dtype=np.uint8)
    rgba[:,:,:3]=np.clip(shade[:,:,None]*255,0,255).astype('uint8')
    rgba[:,:,3]=np.where(finite,255,0).astype('uint8')
    buf=io.BytesIO();Image.fromarray(rgba).save(buf,format='PNG')
    image='data:image/png;base64,'+base64.b64encode(buf.getvalue()).decode('ascii')
    fig=go.Figure()
    if perfil_rio is not None:
        x=perfil_rio.x_utm.to_numpy()-banda['x_B'];y=perfil_rio.y_utm.to_numpy()-banda['y_B']
        valid=(abs(x)<=semilado_m)&(abs(y)<=semilado_m)
        fig.add_trace(go.Scatter(x=np.where(valid,x,np.nan),y=np.where(valid,y,np.nan),mode='lines',
            connectgaps=False,line={'color':'#135479','width':2},name='Río actual',hoverinfo='skip'))
    for mis,color in COLORES_JM.items():
        q=t[t.mis_vis.eq(mis)]
        hover=[texto_terraza(v,incluir_dem=False)+('<br>Dentro del swath 3D' if v.dentro_swath_3D else '<br>Fuera del swath 3D')
               for _,v in q.iterrows()]
        fig.add_trace(go.Scatter(x=q.este_desde_B_m,y=q.norte_desde_B_m,mode='markers+text' if etiquetar else 'markers',
            text=q.terrace_id,textposition='top center',textfont={'color':color,'size':10},
            marker={'color':color,'size':10,'symbol':'diamond' if mis=='5e' else 'circle','line':{'color':'white','width':1}},
            name=f'Terrazas J-M MIS {mis}',legendgroup=f'JM_{mis}',meta={'capa':'terrazas_JM','mis':mis,'n_puntos':len(q)},
            hovertext=hover,hovertemplate='%{hovertext}<extra>J-M 2015</extra>'))
    fig.add_trace(go.Scatter(x=[0],y=[0],mode='markers+text',text=['B inicial'],textposition='bottom center',
        marker={'color':'black','size':10,'symbol':'x'},name='B inicial',hovertemplate='B inicial<extra></extra>'))
    fig.add_shape(type='rect',x0=-banda['ancho_ew_m']/2,x1=banda['ancho_ew_m']/2,
        y0=-banda['largo_ns_m']/2,y1=banda['largo_ns_m']/2,line={'color':'#111111','width':2,'dash':'dash'},fillcolor='rgba(255,255,255,.04)')
    fig.add_layout_image({'source':image,'xref':'x','yref':'y','x':-semilado_m,'y':semilado_m,
        'sizex':2*semilado_m,'sizey':2*semilado_m,'sizing':'stretch','xanchor':'left','yanchor':'top','layer':'below','opacity':.8})
    fig.update_layout(template='plotly_white',height=600,autosize=True,
        title={'text':'Terrazas publicadas cercanas · rectángulo = recorte del swath 3D','font':{'size':16}},
        xaxis={'title':'Oeste (−) / Este (+), m desde B','range':[-semilado_m,semilado_m],'constrain':'domain'},
        yaxis={'title':'Sur (−) / Norte (+), m desde B','range':[-semilado_m,semilado_m],'scaleanchor':'x','scaleratio':1,'constrain':'domain'},
        legend={'orientation':'h','x':.03,'y':-.16},margin={'l':70,'r':35,'t':65,'b':90})
    return fig,t


def exportar_banda(banda, figura, output, caso, abrir=False):
    import webbrowser
    folder=Path(output)
    folder.mkdir(parents=True,exist_ok=True)
    safe=re.sub(r'[^\w.-]+','_',str(caso))
    html=folder/f'{safe}_swath_3D.html'
    script_path=Path(__file__).with_name('seleccionar_puntos_swath.js')
    if not script_path.is_file():
        raise FileNotFoundError(f'Falta el módulo de selección de puntos: {script_path}')
    rid=re.search(r'^R(\d+)',str(caso))
    mis=re.search(r'_(5[ec])_',str(caso))
    selection={'caso_B':str(caso),'rid':int(rid[1]) if rid else None,
               'mis':mis[1] if mis else '', 'x_B':float(banda['x_B']),'y_B':float(banda['y_B']),
               'dem':str(banda['dem'])}
    # Datos de selección explícitos: evita depender del formato binario interno de Plotly.
    selection['east_m']=banda['east_m'].tolist()
    selection['north_m']=banda['north_m'].tolist()
    for key in ['z_m','pendiente_grados']:
        arr=np.asarray(banda[key],float)
        selection[key]=np.where(np.isfinite(arr),arr,None).tolist()
    script=script_path.read_text(encoding='utf8').replace('__CONFIG_SELECCION__',json.dumps(selection,ensure_ascii=False))
    figura.write_html(html,include_plotlyjs=True,full_html=True,
                      config={'scrollZoom':True,'displaylogo':False,'responsive':True},post_script=script)
    for name in ['perfil_NS','perfil_EW']:
        banda[name].to_csv(folder/f'{safe}_{name}.csv',index=False,encoding='utf-8-sig')
    np.savez_compressed(folder/f'{safe}_malla.npz',east_m=banda['east_m'],north_m=banda['north_m'],
                        z_m=banda['z_m'],pendiente_grados=banda['pendiente_grados'])
    meta={k:v for k,v in banda.items() if k not in ['east_m','north_m','z_m','pendiente_grados','perfil_NS','perfil_EW']}
    meta['caso_B']=str(caso)
    (folder/f'{safe}_metadatos.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf8')
    if abrir:
        webbrowser.open(html.resolve().as_uri())
    return html


def indice_html(rows, output):
    from html import escape
    out=Path(output)/'ABRIR_SWATH_3D.html'
    cards=''.join('<li><a href="'+escape(Path(r['html']).name)+'">'+escape(r['caso_B'])+'</a> — '+
                  escape(r['rio'])+' · cobertura '+f"{r['cobertura']:.1%}"+'</li>' for r in rows)
    out.write_text('<!doctype html><html lang="es"><meta charset="utf-8"><title>Bandas 3D alrededor de B</title>'
        '<style>body{font:17px system-ui;max-width:950px;margin:48px auto;line-height:1.7;color:#193347}'
        'a{color:#126aa5}li{margin:14px 0}</style><h1>Relieve y perfiles de banda alrededor de B</h1>'
        '<p>Cada vista contiene relieve 3D y dos perfiles de banda. Azul indica menor pendiente y rojo mayor. '
        'Arrastrar para rotar, rueda para acercar. <strong>Ctrl + clic marca un punto.</strong> '
        'Para guardar cada punto directamente en un archivo, pulsa «Elegir CSV y guardar automáticamente» '
        '(Chrome/Edge). También puedes usar «Descargar CSV». Los archivos funcionan sin internet.</p><ul>'+cards+'</ul>'
        '<p>B sigue siendo tu interpretación inicial. La visualización ayuda a buscar plataforma y paleoacantilado; '
        'no demuestra por sí sola la paleocosta ni modifica K.</p></html>',encoding='utf8')
    return out
