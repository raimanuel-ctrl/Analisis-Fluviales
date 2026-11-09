"""Build the original local review tools with public derived inputs."""
from pathlib import Path
import base64
import gzip
import json
import sys
from datetime import datetime, timezone

BASE=Path(__file__).resolve().parents[1]


def generar_explorador(output=None):
    source=BASE/'herramientas/explorador_kp'
    raw=(source/'datos_visores.json').read_bytes()
    data=json.loads(raw)
    assert data['initial']['mediciones']==[]
    payload=base64.b64encode(gzip.compress(raw,compresslevel=7,mtime=0)).decode('ascii')
    panels=''
    for i,label in enumerate(['Panel A','Panel B']):
        panels+=f'''<section class="panel"><div class="panelhead"><strong>{label}</strong><select id="p{i}Type" aria-label="Tipo de perfil {label}"><option value="long">Longitudinal</option><option value="chi">χ–elevación</option><option value="tau">τ–elevación</option></select><div><button id="p{i}Reset">Perfil completo</button><button id="p{i}PNG">PNG</button></div></div><div class="regctrl"><label><input id="p{i}Reg" type="checkbox">Regresión en χ</label><select id="p{i}Scope" aria-label="Alcance de ajuste {label}"><option value="full">Perfil completo</option><option value="3000">3 km por lado</option></select><span>Sin zoom</span><label>Excluir ±<input id="p{i}Gap" type="number" value="100" min="0" max="5000" step="25" aria-label="Exclusion en metros {label}">m</label><button id="p{i}Cross" disabled>Guardar cruce</button></div><div class="canvasbox" id="p{i}Box"><canvas id="p{i}Base"></canvas><canvas class="overlay" id="p{i}Overlay" tabindex="0" aria-label="Perfil interactivo {label}"></canvas></div><div class="regstats" id="p{i}Stats">Mueve el cursor cerca del cauce.</div></section>'''
    text=(source/'index.template.html').read_text(encoding='utf-8')
    for key,value in [('__PANELS__',panels),('__MODEL__',(source/'model.js').read_text(encoding='utf-8')),
                      ('__PAYLOAD__',payload),('__UI__',(source/'ui.js').read_text(encoding='utf-8'))]:
        text=text.replace(key,value)
    dest=Path(output) if output else BASE/'visores/Explorador_KP.html'
    dest.parent.mkdir(exist_ok=True,parents=True);dest.write_text(text,encoding='utf-8')
    return dest


def generar_b_3d(output=None, rids=None, paso_m=None, semiancho_ew_m=1000.,semilargo_ns_m=2000.):
    """No original spatial terrace catalogue is included in the public viewer.

    Read a compatible, authorized local CSV by configuring terrazas_opcional.
    An HTML regenerated with that CSV embeds its points and must remain local.
    """
    import pandas as pd
    from . import analisis as a
    support=BASE/'herramientas/visor_b_3d'
    if str(support) not in sys.path:sys.path.insert(0,str(support))
    import swath_b_3d as sb
    import visor_rios_swath as vr
    cfg=json.loads((BASE/'configuracion.json').read_text(encoding='utf-8'))
    dem=Path(cfg['dem'])
    if not dem.is_absolute():dem=BASE/dem
    if not dem.is_file() or dem.stat().st_size<1000:
        raise FileNotFoundError('Descarga el DEM mediante Git LFS y revisa configuracion.json antes de regenerar el relieve.')
    cols=['terrace_id','este_utm18s_m','norte_utm18s_m','mis_asignado','elevacion_m',
          'error_elevacion_m','edad_referencia_ka','incertidumbre_edad_ka']
    pts=pd.DataFrame(columns=cols)
    if cfg.get('terrazas_opcional'):
        t=Path(cfg['terrazas_opcional']);t=t if t.is_absolute() else BASE/t
        pts=pd.read_csv(t)
    cases=[]
    for b in a.leer('selecciones_b.csv').itertuples():
        if rids is not None and b.rid not in rids:continue
        band=sb.preparar_banda(dem,b.x_b_m,b.y_b_m,semiancho_ew_m,semilargo_ns_m,paso_m=paso_m)
        profile=a.perfil(b.rid)
        fig=sb.graficar_banda(band,b.rio,b.id_b,perfil_rio=profile,
            puntos_JM=pts if not pts.empty else None,exageracion_vertical=1.)
        context,context_pts=sb.mapa_contexto_terrazas(band,pts,perfil_rio=profile,semilado_m=3000.)
        if pts.empty:context.update_layout(title='Contexto del relieve y cauce · rectángulo = recorte del swath 3D')
        cases.append(dict(banda=band,figura=fig,caso_B=b.id_b,rid=b.rid,rio=b.rio,mis=b.mis,
                          contexto=context,terrazas_contexto=context_pts,perfil_rio=profile,
                          semilado_contexto_m=3000.,exageracion_vertical=1.))
    if not cases:raise ValueError('No hay casos B para los rios solicitados.')
    dest=Path(output) if output else BASE/'resultados'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')/'Visor_B_3D.html'
    dest.parent.mkdir(exist_ok=True,parents=True)
    out=vr.exportar_visor(cases,dest.parent,nombre=dest.name,mostrar_5e=not pts.empty,mostrar_5c=not pts.empty)
    if pts.empty:
        text=out.read_text(encoding='utf-8').replace('<header><h1>','<header><p><b>Versión pública:</b> relieve, cauces y B adoptados. El catálogo espacial de terrazas se aporta externamente. Los clics registran candidatos de revisión.</p><h1>')
        out.write_text(text,encoding='utf-8')
    return out


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('visor',choices=['kp','b3d'])
    p.add_argument('--output',type=Path)
    p.add_argument('--rio',type=int,action='append')
    args=p.parse_args()
    print(generar_explorador(args.output) if args.visor=='kp' else generar_b_3d(args.output,rids=args.rio))
