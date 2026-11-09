"""HTML único con selectores río/escenario y puntos independientes por caso."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import json
import re
import base64
import gzip
import numpy as np
from plotly.offline import get_plotlyjs


def figura_json(fig):
    # Arrays explícitos: también editables por el controlador del navegador.
    def plain(v):
        if isinstance(v,dict) and 'bdata' in v and 'dtype' in v:
            a=np.frombuffer(base64.b64decode(v['bdata']),dtype=v['dtype'])
            if v.get('shape'):a=a.reshape(tuple(int(n) for n in v['shape'].split(',')))
            return np.where(np.isfinite(a),a,None).tolist()
        if isinstance(v,dict):return {k:plain(x) for k,x in v.items()}
        if isinstance(v,list):return [plain(x) for x in v]
        return v
    return plain(json.loads(fig.to_json()))


def raster_para_ajuste(banda, semilado_m=3000.):
    """DEM nativo y pendiente con halo. Float32/base64, lectura local en el HTML."""
    import base64
    import rasterio
    from rasterio.windows import from_bounds
    x,y=banda['x_B'],banda['y_B']
    with rasterio.open(banda['dem']) as src:
        halo=3*max(src.res)
        win=from_bounds(x-semilado_m-halo,y-semilado_m-halo,x+semilado_m+halo,y+semilado_m+halo,
                        src.transform).round_offsets().round_lengths()
        z=src.read(1,window=win,masked=True,boundless=True).astype(float).filled(np.nan)
        gy,gx=np.gradient(z,abs(src.res[1]),abs(src.res[0]))
        slope=np.degrees(np.arctan(np.hypot(gx,gy)))
        tr=src.window_transform(win)
        def pack(a):return base64.b64encode(np.asarray(a,dtype='<f4').tobytes(order='C')).decode('ascii')
        return {'z_b64':pack(z),'pendiente_b64':pack(slope),'nx':z.shape[1],'ny':z.shape[0],
                'x0':tr.c+tr.a/2,'y0':tr.f+tr.e/2,'dx':tr.a,'dy':tr.e,
                'resolucion_m':max(abs(src.res[0]),abs(src.res[1])),
                'max_ancho_m':2*semilado_m,'max_largo_m':2*semilado_m}


def configuracion(banda, caso, rid, rio, mis):
    cfg={'caso_B':str(caso),'rid':int(rid),'rio':str(rio),'mis':str(mis),
         'x_B':float(banda['x_B']),'y_B':float(banda['y_B']),'dem':str(banda['dem']),
         'east_m':banda['east_m'].tolist(),'north_m':banda['north_m'].tolist()}
    for key in ['z_m','pendiente_grados']:
        arr=np.asarray(banda[key],float)
        cfg[key]=np.where(np.isfinite(arr),arr,None).tolist()
    return cfg


def exportar_visor(casos, output, nombre='Visor_rios_swath_3D.html', mostrar_5e=True, mostrar_5c=True):
    """casos: lista de banda, figura, caso_B, rid, rio, mis. Todo queda dentro de un HTML."""
    if not casos:
        raise ValueError('No hay bandas para el visor.')
    ids=[c['caso_B'] for c in casos]
    if len(ids)!=len(set(ids)):
        raise ValueError('caso_B repetido en el visor.')
    payload=[]
    counts={}
    for item in casos:
        key=(int(item['rid']),str(item['mis']))
        counts[key]=counts.get(key,0)+1
        cfg=configuracion(item['banda'],item['caso_B'],item['rid'],item['rio'],item['mis'])
        cfg['JM_en_3D']={mis:sum(int(t.meta.get('n_puntos',0)) for t in item['figura'].data
                                if isinstance(t.meta,dict) and t.meta.get('capa')=='terrazas_JM' and t.meta.get('mis')==mis)
                         for mis in ['5e','5c']}
        context=item.get('contexto')
        cfg['JM_en_contexto']={mis:sum(int(t.meta.get('n_puntos',0)) for t in context.data
                                      if isinstance(t.meta,dict) and t.meta.get('capa')=='terrazas_JM' and t.meta.get('mis')==mis)
                               for mis in ['5e','5c']} if context is not None else None
        pts=item.get('terrazas_contexto')
        cfg['JM_catalogo']=(pts[['terrace_id','mis_vis','este_utm18s_m','norte_utm18s_m','elevacion_m',
            'error_elevacion_m','edad_referencia_ka','incertidumbre_edad_ka']].to_dict('records') if pts is not None else [])
        p=item.get('perfil_rio')
        if p is not None:
            mask=(abs(p.x_utm-cfg['x_B'])<=item.get('semilado_contexto_m',3000)) & (abs(p.y_utm-cfg['y_B'])<=item.get('semilado_contexto_m',3000))
            cfg['river']={k:np.where(mask,p[col].to_numpy()-offset,None).tolist() for k,col,offset in
                [('x','x_utm',cfg['x_B']),('y','y_utm',cfg['y_B']),('z','elevacion_m',0)]}
        cfg['ancho_ew_m']=float(item['banda']['ancho_ew_m'])
        cfg['largo_ns_m']=float(item['banda']['largo_ns_m'])
        cfg['exageracion_vertical']=float(item.get('exageracion_vertical',1))
        # La malla se reconstruye desde la reserva nativa al abrir/aplicar vista.
        # Evitar guardar varias copias de los mismos cientos de miles de nodos.
        for name in ['east_m','north_m','z_m','pendiente_grados']:cfg.pop(name,None)
        figjson=figura_json(item['figura'])
        for tr in figjson['data']:
            if tr.get('meta',{}).get('capa'):
                for name in ['x','y','z','surfacecolor','customdata','text','hovertext']:
                    if name in tr:tr[name]=[]
        payload.append({'cfg':cfg,'etiqueta':f"MIS {item['mis']} · B inicial {counts[key]}",
                        'fig':figjson,
                        'contexto':figura_json(context) if context is not None else None,
                        'raster':raster_para_ajuste(item['banda'],item.get('semilado_contexto_m',3000))})
    # La configuración inicial efectiva puede recuperarse del último caso visto.
    # Instalar con el contexto resuelto por el controlador, no uno fijado en Python.
    selector=Path(__file__).with_name('seleccionar_puntos_swath.js').read_text(encoding='utf8')
    selector=selector.replace("document.getElementById('{plot_id}')","document.getElementById('swath-map')")
    selector=selector.replace('__CONFIG_SELECCION__','initial.cfg')
    controller=Path(__file__).with_name('controlador_visor_rios.js').read_text(encoding='utf8')
    controller=controller.replace('__INSTALAR_SELECTOR__',selector)
    dynamic=Path(__file__).with_name('ajustar_vista_swath.js').read_text(encoding='utf8')
    controller=dynamic+'\n'+controller
    data=json.dumps({'casos':payload,'mostrarTerrazas':{'5e':bool(mostrar_5e),'5c':bool(mostrar_5c)},
                    'plotConfig':{'scrollZoom':True,'displaylogo':False,'responsive':True}},
                    ensure_ascii=False,separators=(',',':'),allow_nan=False).replace('<',chr(92)+'u003c')
    data=base64.b64encode(gzip.compress(data.encode('utf8'),compresslevel=6,mtime=0)).decode('ascii')
    nrios=len(set(int(c['rid']) for c in casos))
    html='''<!doctype html><html lang="es"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Ríos · swath 3D y selección de puntos</title>
<style>body{margin:0;color:#19374b;background:white;font:14px system-ui,sans-serif}
header{margin:12px;padding:14px 18px;background:#e9f1f7;border:1px solid #bfd1df;border-radius:8px}
h1{font-size:21px;margin:0 0 9px}.controles{display:flex;align-items:center;gap:10px;flex-wrap:wrap}
label{display:flex;gap:7px;align-items:center;font-weight:600}select,button{font:inherit;padding:7px 10px;border:1px solid #9bb5c8;border-radius:5px;background:white;color:#164b70}
button{cursor:pointer}button.principal{background:#176c9e;color:white}button:disabled,select:disabled{opacity:.6}
.ayuda{font-size:12px;color:#506a7c;margin-top:9px;line-height:1.5}#visor-status{margin-top:7px;font-size:13px}
#visor-total{font-weight:600;font-size:13px}#swath-map{width:100%;min-height:850px}</style></head><body>
<header><h1>Relieve 3D y puntos por río</h1><div class="controles">
<label>Río <select id="visor-rio" aria-label="Río" disabled></select></label>
<label>MIS / escenario B <select id="visor-escenario" aria-label="MIS y escenario B" disabled></select></label>
<button id="visor-vincular-global" class="principal">Elegir CSV de todos los ríos</button>
<button id="visor-descargar-todos">Descargar todos (CSV)</button>
<button id="visor-importar">Importar puntos CSV</button><input id="visor-archivo-csv" type="file" accept=".csv,text/csv" multiple hidden>
<span id="visor-total"></span></div>
<div class="controles" style="margin-top:9px">
<label style="color:#c10caf"><input type="checkbox" id="visor-mis5e" aria-label="Mostrar terrazas MIS 5e" checked>Terrazas J-M MIS 5e</label>
<label style="color:#138c45"><input type="checkbox" id="visor-mis5c" aria-label="Mostrar terrazas MIS 5c" checked>Terrazas J-M MIS 5c</label>
<span id="visor-terrazas-info" style="font-size:12px"></span></div>
<div class="controles" style="margin-top:10px">
<label>Ancho O–E (m) <input type="number" id="visor-ancho" aria-label="Ancho oeste este en metros" value="1000" min="100" max="6000" step="50" style="width:85px"></label>
<label>Largo N–S (m) <input type="number" id="visor-largo" aria-label="Largo norte sur en metros" value="2000" min="100" max="6000" step="50" style="width:85px"></label>
<label>Exageración vertical <input type="number" id="visor-ev" aria-label="Exageración vertical" value="1" min="0.1" max="30" step="0.5" style="width:70px"></label>
<label><input type="checkbox" id="visor-color-pendiente" checked> Colorear pendiente</label><button id="visor-aplicar">Aplicar vista</button><span id="visor-malla-info" style="font-size:12px"></span></div>
<div id="visor-status" role="status" aria-live="polite">Preparando las vistas…</div>
<div class="ayuda">__NRIOS__ ríos con B manual · __NCASOS__ escenarios. Cambia de río aquí; cada escenario conserva sus puntos.
Ctrl + clic sobre el relieve para marcar. El CSV de todos incluye río, MIS y escenario; el panel inferior permite guardar sólo el escenario actual.
Usa «Importar puntos CSV» para recuperar archivos anteriores o trasladar selecciones a otro navegador.
El respaldo local pertenece al navegador y origen usados: conserva también un CSV. Azul = pendiente suave; rojo = fuerte. B y K no se modifican automáticamente.
Los símbolos J-M se apoyan en el DEM con un despegue visual de 3 m; su altura publicada y edad MIS aparecen al señalar el punto.
El mapa inferior muestra también terrazas cercanas que quedan fuera del swath 3D.</div></header>
<div id="swath-map"></div><section id="visor-contexto-wrap"><div id="visor-contexto" style="width:100%;min-height:600px"></div></section><script>__PLOTLY__</script>
<script type="application/json" id="datos-visor">__DATOS__</script><script>__CONTROLADOR__</script></body></html>'''
    html=html.replace('__NRIOS__',str(nrios)).replace('__NCASOS__',str(len(casos)))
    html=html.replace('__PLOTLY__',get_plotlyjs()).replace('__DATOS__',data).replace('__CONTROLADOR__',controller)
    target=Path(output)/nombre
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(html,encoding='utf8')
    return target
