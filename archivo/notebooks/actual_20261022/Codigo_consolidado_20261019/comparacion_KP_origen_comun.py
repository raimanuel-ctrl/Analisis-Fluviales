"""Comparación de elevaciones de KP y tiempos condicionales de propagación.

No recalibra K ni escoge escenarios. La elevación nunca se desplaza.
"""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from shapely.geometry import LineString, Point


def _bool(s):
    return s.astype(str).str.lower().isin(['true', '1'])


def _rango(s, cifras=1):
    a = pd.to_numeric(pd.Series(s), errors='coerce').dropna()
    if a.empty:
        return '—'
    lo, hi = float(a.min()), float(a.max())
    return f'{lo:.{cifras}f}' if np.isclose(lo, hi) else f'{lo:.{cifras}f}–{hi:.{cifras}f}'


def _color(rid):
    order = list(range(0, 20, 2)) + list(range(1, 20, 2))
    return plt.get_cmap('tab20')(order[(int(rid)-1) % 20])


def _nombre(s):
    return str(s).replace('Estero ', 'E. ').replace('Rio ', 'R. ')


def _etiquetas(ax, puntos):
    """Coloca rótulos cortos sin solaparlos, con líneas al KP real."""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    occupied = []
    limits = ax.get_window_extent(renderer)
    offsets = [(8, 9), (8, -18), (-35, 10), (-35, -18), (15, 25),
               (-40, 28), (20, -33), (-42, -35), (10, 43), (-40, 47),
               (48, 6), (-70, 6), (35, 55), (-55, 64), (30, -57)]
    for x, y, text, color in sorted(puntos, key=lambda v: (v[0], v[1])):
        best = None
        for ox, oy in offsets:
            note = ax.annotate(text, (x, y), xytext=(ox, oy), textcoords='offset points',
                               fontsize=9, weight='bold', color='#162e3c', zorder=20,
                               bbox=dict(boxstyle='round,pad=.20', fc='white', ec=color, alpha=.94),
                               arrowprops=dict(arrowstyle='-', color=color, lw=.65))
            box = note.get_bbox_patch().get_window_extent(renderer).expanded(1.10, 1.15)
            score = sum(box.overlaps(b) for b in occupied)
            score += 50 * (not limits.contains(box.x0, box.y0) or not limits.contains(box.x1, box.y1))
            if best is None or score < best[0]:
                best = (score, ox, oy, box)
            note.remove()
            if score == 0:
                break
        _, ox, oy, box = best
        ax.annotate(text, (x, y), xytext=(ox, oy), textcoords='offset points',
                    fontsize=9, weight='bold', color='#162e3c', zorder=20,
                    bbox=dict(boxstyle='round,pad=.20', fc='white', ec=color, alpha=.94),
                    arrowprops=dict(arrowstyle='-', color=color, lw=.65))
        occupied.append(box)


def _tabla(ax, rows, headers, widths, rids, fontsize=10, bbox=(0, .08, 1, .82)):
    ax.axis('off')
    t = ax.table(cellText=rows, colLabels=headers, colWidths=widths,
                 cellLoc='center', loc='center', bbox=bbox)
    t.auto_set_font_size(False)
    t.set_fontsize(fontsize)
    for (i, j), cell in t.get_celld().items():
        cell.set_edgecolor('#c4d0db')
        cell.set_linewidth(.5)
        if i == 0:
            cell.set_facecolor('#203d51'); cell.set_text_props(color='white', weight='bold')
        else:
            cell.set_facecolor('#f0f5f9' if i % 2 else 'white')
            if j == 0:
                cell.set_facecolor(_color(rids[i-1]))
                cell.set_text_props(color='#102432', weight='bold')
    return t


def generar(base, salida, A0=1e6, theta=.40, rios=None,
            mostrar_todas=True, mostrar_rango=True, limite_tau_ka=800.,
            limite_altura_m=None, interactivo=True):
    """Lee sólo una corrida consolidada. rios=None incluye los 18 cauces.

    χ se pone a cero en el outlet analítico de cada río (excluyendo extensiones
    recuperadas para reconstruir B). τ usa cada B, sin desplazar la elevación.
    Los límites sólo recortan la vista; las tablas conservan todos los valores.
    """
    base, salida = Path(base), Path(salida)
    out = salida/'Comparacion_KP_origen_comun'
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((base/'Configuracion.json').read_text(encoding='utf-8-sig'))
    if not np.isclose(theta, cfg['theta']) or not np.isclose(A0, cfg['A0_m2']):
        raise ValueError('A0 y theta deben coincidir con la corrida de K; no cambiar sólo el gráfico.')
    k = pd.read_csv(salida/'K_todos_los_cruces.csv', float_precision='round_trip')
    kt = pd.read_csv(salida/'Tau/Tau_KP_por_escenario.csv', float_precision='round_trip')
    manual_path = base/cfg['kp_manual_csv']
    manual = pd.read_csv(manual_path, sep=cfg.get('kp_separador', ';'), decimal=cfg.get('kp_decimal', ','))
    nombres = pd.read_csv(base/'Productos/03_K_U/v2/resultados_K_U_v1.csv').set_index('rid').rio.to_dict()
    ids = sorted(nombres if rios is None else set(map(int, rios)))
    k = k[k.rid.isin(ids) & np.isfinite(k.K_m02_anio) & k.K_m02_anio.gt(0)].copy()
    kt = kt[kt.rid.isin(ids)].copy()
    perfiles, puntos, origenes = {}, [], []
    for rid in ids:
        p = pd.read_csv(salida/f'perfil_R{rid:02d}_utilizado.csv', float_precision='round_trip')
        if not np.all(np.diff(p.chi_m) >= 0):
            raise ValueError(f'χ no monótono: R{rid:02d}')
        original = p[~_bool(p.tramo_recuperado_aguas_abajo)] if 'tramo_recuperado_aguas_abajo' in p else p
        if original.empty:
            raise ValueError(f'Sin outlet analítico original: R{rid:02d}')
        o = original.iloc[0]
        p = original.copy()
        p['chi_rel_m'] = p.chi_m-float(o.chi_m)
        assert p.chi_rel_m.iloc[0] == 0
        perfiles[rid] = p
        origenes.append(dict(rid=rid, rio=nombres[rid], chi_outlet_analisis_m=float(o.chi_m),
                             z_outlet_analisis_m=float(o.elevacion_m), x_utm=float(o.x_utm), y_utm=float(o.y_utm)))
        line = LineString(p[['x_utm', 'y_utm']].to_numpy())
        sgeom = np.r_[0, np.cumsum(np.hypot(np.diff(p.x_utm), np.diff(p.y_utm)))]
        rows = manual[pd.to_numeric(manual.id, errors='coerce').eq(rid)].dropna(subset=['id_kp', 'x_utm', 'y_utm'])
        for num, (_, v) in enumerate(rows.iterrows(), 1):
            point = Point(float(v.x_utm), float(v.y_utm)); s = line.project(point)
            distance = float(point.distance(line))
            chi = float(np.interp(s, sgeom, p.chi_m))
            z = float(np.interp(s, sgeom, p.elevacion_m))
            label = f'R{rid:02d}' + (f'-KP{num}' if len(rows)>1 else '')
            puntos.append(dict(rid=rid, rio=nombres[rid], rotulo=label, id_kp=str(v.id_kp),
                               chi_KP_m=chi, chi_rel_KP_m=chi-float(o.chi_m), z_KP_m=z,
                               separacion_cauce_m=distance, valido=distance<=35))
    kp = pd.DataFrame(puntos)
    # Las edades se vuelven a comprobar con la fórmula, sin sustituir valores.
    usable = kt[kt.estado.eq('utilizable') & np.isfinite(kt.tau_KP_ka)].copy()
    expected = (usable.chi_kp_m-usable.chi_B_m)/(usable.K_m02_anio*A0**theta)/1000
    assert np.allclose(expected, usable.tau_KP_ka, rtol=1e-12, atol=1e-10)
    differences = []
    for _, v in usable.iterrows():
        match = kp[kp.id_kp.eq(v.id_kp) & kp.rid.eq(v.rid)]
        if len(match) != 1 or not match.iloc[0].valido:
            raise ValueError(f'KP de τ no reproducible: {v.id_kp}')
        differences.append(abs(float(match.iloc[0].z_KP_m)-v.z_kp_m))
        assert np.isclose(float(match.iloc[0].chi_KP_m), v.chi_kp_m, atol=1e-7, rtol=0)
    assert not differences or max(differences)<1e-7
    resumen = []
    for rid in ids:
        m = kp[kp.rid.eq(rid) & kp.valido]
        # Se mantiene cada KP si en el futuro se selecciona más de uno por río.
        for _, mark in (m.iterrows() if len(m) else [(None, None)]):
            q = usable[usable.rid.eq(rid)]
            if mark is not None:
                q = q[q.id_kp.eq(mark.id_kp)]
            resumen.append(dict(rid=rid, rio=nombres[rid], rotulo=mark.rotulo if mark is not None else f'R{rid:02d}',
                 id_kp=mark.id_kp if mark is not None else None, z_KP_m=mark.z_KP_m if mark is not None else np.nan,
                 tau_5e_min_ka=q.loc[q.mis.eq('5e'), 'tau_KP_ka'].min(), tau_5e_max_ka=q.loc[q.mis.eq('5e'), 'tau_KP_ka'].max(),
                 tau_5c_min_ka=q.loc[q.mis.eq('5c'), 'tau_KP_ka'].min(), tau_5c_max_ka=q.loc[q.mis.eq('5c'), 'tau_KP_ka'].max(),
                 n_escenarios_tau=len(q), estado='Sin KP manual' if mark is None else ('τ calculado' if len(q) else 'Sin K calibrado')))
    resumen = pd.DataFrame(resumen)
    por_B = []
    for (rid, case), q in k.groupby(['rid', 'caso_B'], sort=True):
        t = usable[usable.caso_B.eq(case)]
        for kid in (t.id_kp.unique().tolist() if len(t) else [None]):
            tt = t[t.id_kp.eq(kid)] if kid is not None else t
            por_B.append(dict(rid=rid, rio=q.rio.iloc[0], caso_B=case, id_B=q.id_B_corto.iloc[0], mis=q.mis.iloc[0],
                id_kp=kid, n_K=len(q), n_K_con_tau=len(tt), K_min_m02_anio=q.K_m02_anio.min(), K_max_m02_anio=q.K_m02_anio.max(),
                tau_KP_min_ka=tt.tau_KP_ka.min(), tau_KP_max_ka=tt.tau_KP_ka.max(),
                z_KP_m=tt.z_kp_m.iloc[0] if len(tt) else np.nan,
                estado='τ calculado' if len(tt) else 'Sin KP utilizable'))
    por_B = pd.DataFrame(por_B)
    for name, data in [('Resumen_KP_por_rio', resumen), ('Resumen_KP_por_B', por_B),
                       ('KP_en_perfiles_actuales', kp), ('Origenes_analiticos_chi', pd.DataFrame(origenes))]:
        data.to_csv(out/f'{name}.csv', index=False, encoding='utf-8-sig')
    fig = plt.figure(figsize=(23, 11), facecolor='white')
    gs = fig.add_gridspec(1, 2, width_ratios=[1.38, 1.], left=.055, right=.985, top=.88, bottom=.14, wspace=.09)
    ax = fig.add_subplot(gs[0]); table_ax = fig.add_subplot(gs[1])
    ymax = max(float(p.elevacion_m.max()) for p in perfiles.values())*1.06
    for rid, p in perfiles.items():
        ax.plot(p.chi_rel_m, p.elevacion_m, color=_color(rid), lw=1.35, alpha=.85)
    valid = kp[kp.valido]
    ax.scatter(valid.chi_rel_KP_m, valid.z_KP_m, c=[_color(i) for i in valid.rid], s=155, marker='*', edgecolors='black', linewidths=.7, zorder=9)
    ax.set(xlabel='χ − χ(outlet de análisis) (m) · aguas arriba →', ylabel='Elevación DEM (m) · sin desplazamiento vertical',
           xlim=(0, max(p.chi_rel_m.max() for p in perfiles.values())*1.02), ylim=(0, limite_altura_m or ymax))
    ax.set_title('Perfiles χ superpuestos · origen horizontal alineado', loc='left', fontsize=15, pad=15)
    ax.grid(alpha=.20)
    ax.ticklabel_format(axis='x', style='plain', useOffset=False)
    _etiquetas(ax, [(v.chi_rel_KP_m, v.z_KP_m, v.rotulo, _color(v.rid)) for _, v in valid.iterrows() if v.z_KP_m<=ax.get_ylim()[1]])
    rows = []
    for _, v in resumen.iterrows():
        tf = lambda mis: _rango([v[f'tau_{mis}_min_ka'], v[f'tau_{mis}_max_ka']]) if np.isfinite(v[f'tau_{mis}_min_ka']) else ('Sin KP' if pd.isna(v.id_kp) else 'Sin K')
        rows.append([v.rotulo, _nombre(v.rio), f'{v.z_KP_m:.1f}' if np.isfinite(v.z_KP_m) else '—', tf('5e'), tf('5c')])
    _tabla(table_ax, rows, ['KP / río', 'Nombre', 'Altura KP\n(m)', 'τ con MIS 5e\n(ka)', 'τ con MIS 5c\n(ka)'],
           [.12, .32, .13, .215, .215], resumen.rid.tolist(), fontsize=9.8, bbox=(0,.04,1,.88))
    table_ax.set_title('Tiempos de propagación en los KP manuales', fontsize=14, loc='left', pad=15)
    fig.suptitle('Alturas de los knickpoints y tiempos estimados con K propios', fontsize=21, x=.055, ha='left', y=.975, weight='bold', color='#203d51')
    fig.text(.055,.924,'χ = 0 en el inicio de cada cauce analizado; los orígenes no coinciden geográficamente. Estrellas = KP manuales.', fontsize=12)
    fig.text(.055,.075,'Tabla: τKP = [χ(KP) − χ(B)] / (K · A₀^0,40) / 1000. Intervalos = mínimo–máximo de los escenarios B/K de cada MIS.', fontsize=12)
    fig.text(.055,.042,'MIS identifica la terraza que calibra K, no la edad del KP. τ es un tiempo modelado, con K uniforme y área actual; no una datación independiente.', fontsize=12)
    fig.savefig(out/'Perfiles_chi_origen_comun_KP.png', dpi=180, facecolor='white')
    plt.close(fig)
    # Segunda vista: todas las calibraciones válidas en una coordenada temporal.
    fig = plt.figure(figsize=(21, 10), facecolor='white')
    gs = fig.add_gridspec(1, 2, width_ratios=[1.42,1.], left=.055,right=.985,top=.86,bottom=.20,wspace=.10)
    ax = fig.add_subplot(gs[0]); ta = fig.add_subplot(gs[1])
    tau_puntos = []
    curves_for_html = []
    for (rid, case), q in k.groupby(['rid','caso_B'], sort=True):
        pfull = pd.read_csv(salida/f'perfil_R{int(rid):02d}_utilizado.csv', float_precision='round_trip')
        b = q.iloc[0]; tail = pfull[pfull.chi_m>b.chi_B_m+1e-8]
        chi = np.r_[b.chi_B_m,tail.chi_m.to_numpy()]; z=np.r_[b.z_B_m,tail.elevacion_m.to_numpy()]
        ts = np.array([(chi-v.chi_B_m)/(v.K_m02_anio*A0**theta)/1000 for _,v in q.iterrows()])
        assert np.all(ts[:,0] == 0) and np.all(ts>=0)
        style = '-' if b.mis=='5e' else '--'
        if mostrar_rango and len(ts)>1:
            ax.fill_betweenx(z, ts.min(axis=0), ts.max(axis=0), color=_color(rid), alpha=.08)
        for (_,v),tau in zip(q.iterrows(),ts):
            if mostrar_todas:
                ax.plot(tau,z,color=_color(rid),ls=style,lw=.75,alpha=.24)
            curves_for_html.append((int(rid), v, tau, z))
        if not mostrar_todas or mostrar_rango:
            ax.plot(ts.min(axis=0),z,color=_color(rid),ls=style,lw=.85,alpha=.6)
            ax.plot(ts.max(axis=0),z,color=_color(rid),ls=style,lw=.85,alpha=.6)
    temporal_rows=[]; temporal_rids=[]
    for (rid, mis), q in k.groupby(['rid','mis'], sort=True):
        tt = usable[usable.rid.eq(rid)&usable.mis.eq(mis)]
        shape='*' if mis=='5e' else 'D'
        if len(tt):
            lo,hi=tt.tau_KP_ka.min(),tt.tau_KP_ka.max();zz=tt.z_kp_m.iloc[0]
            ax.hlines(zz,lo,hi,color=_color(rid),lw=3,zorder=8)
            ax.scatter(tt.tau_KP_ka,tt.z_kp_m,marker=shape,s=60,c=[_color(rid)],edgecolor='#1c2930',lw=.45,zorder=9)
            tau_puntos.append(((lo+hi)/2,zz,f'R{int(rid):02d} · {mis}',_color(rid)))
        temporal_rows.append([f'R{int(rid):02d}', mis,str(len(q)),
             _rango(q.K_m02_anio*1e6,2),f'{tt.z_kp_m.iloc[0]:.1f}' if len(tt) else '—',_rango(tt.tau_KP_ka) if len(tt) else 'Sin KP'])
        temporal_rids.append(rid)
    ax.set(xlabel='τ relativo a cada B (ka) · aguas arriba →',ylabel='Elevación DEM (m) · sin desplazamiento vertical',
           xlim=(0, limite_tau_ka if limite_tau_ka is not None else max(v[2][-1] for v in curves_for_html)*1.02),
           ylim=(0, limite_altura_m or ymax))
    ax.grid(alpha=.20)
    ax.set_title('Perfiles τ superpuestos · cada B se alinea en τ = 0', loc='left',fontsize=14,pad=15)
    _etiquetas(ax, [v for v in tau_puntos if v[0]<=ax.get_xlim()[1] and v[1]<=ax.get_ylim()[1]])
    handles=[Line2D([0],[0],color=_color(r),label=f'R{int(r):02d} · {_nombre(nombres[r])}') for r in sorted(k.rid.unique())]
    handles += [Line2D([0],[0],color='#333333',ls='-',marker='*',label='MIS 5e'),Line2D([0],[0],color='#333333',ls='--',marker='D',label='MIS 5c')]
    ax.legend(handles=handles,fontsize=9,loc='upper left',ncol=2,framealpha=.94)
    _tabla(ta,temporal_rows,['Río','MIS','N K','K × 10⁶\n(m⁰·²/año)','Altura KP\n(m)','τ en KP\n(ka)'],
           [.12,.08,.09,.25,.15,.31],temporal_rids,fontsize=10.2,bbox=(0,.31,1,.61))
    ta.set_title('Rangos de escenarios; no intervalos de confianza',loc='left',fontsize=13,pad=15)
    ta.text(.01,.24,'N K = número de calibraciones B–P_T conservadas.\nAltura = elevación del KP en el perfil actual.\nCada fila agrupa los distintos B y P_T del mismo MIS.\nChanco conserva sus K, pero no tiene KP manual.',transform=ta.transAxes,fontsize=11,va='top',linespacing=1.6)
    fig.suptitle('Comparación temporal con todos los K calibrados',fontsize=21,x=.055,ha='left',y=.97,weight='bold',color='#203d51')
    fig.text(.055,.915,'Color = río; línea continua/estrella = MIS 5e; línea discontinua/rombo = MIS 5c. Los KP conservan su elevación.',fontsize=12)
    zoom = f'La vista termina en {limite_tau_ka:g} ka; el cálculo y las curvas interactivas conservan el perfil completo.' if limite_tau_ka is not None else 'La vista incluye el perfil completo.'
    fig.text(.055,.13,zoom,fontsize=12)
    fig.text(.055,.09,'τ no es una datación independiente del KP: extrapola cada K de B–P_T a la posición del KP; τ(P_T) = T por construcción.',fontsize=12)
    fig.text(.055,.05,'Las posiciones de B y los diferentes K producen un abanico de escenarios; no se ha elegido un K representativo ni una edad final.',fontsize=12)
    fig.savefig(out/'Perfiles_tau_origen_B_KP.png',dpi=180,facecolor='white');plt.close(fig)
    html = None
    if interactivo:
        html = _interactivo(out, perfiles, kp, resumen, curves_for_html, usable, nombres, limite_tau_ka, limite_altura_m, mostrar_todas)
    verify={'n_perfiles':len(perfiles),'n_KP_manuales_proyectados':len(kp),'n_KP_validos':int(kp.valido.sum()),
            'n_K':len(k),'n_K_con_tau_KP':len(usable),'n_B':int(k.caso_B.nunique()),
            'max_error_z_KP_vs_corrida_tau_m':max(differences) if differences else None,
            'formula_tau':'(chi_KP-chi_B)/(K*A0**theta)/1000 [ka]',
            'origen_chi':'outlet analítico individual; se excluye la extensión recuperada para B',
            'origen_tau':'cada B individual; no se desplazan alturas',
            'alcance':'tiempo condicional de propagación; no edad tectónica independiente',
            'A0_m2':A0,'theta':theta,'limite_tau_solo_visual_ka':limite_tau_ka,
            'fuentes_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [salida/'K_todos_los_cruces.csv',salida/'Tau/Tau_KP_por_escenario.csv',manual_path]}}
    (out/'Verificacion_comparacion_KP.json').write_text(json.dumps(verify,ensure_ascii=False,indent=2),encoding='utf-8')
    return {'png_chi':out/'Perfiles_chi_origen_comun_KP.png','png_tau':out/'Perfiles_tau_origen_B_KP.png',
            'html':html,'resumen':resumen,'por_B':por_B,'verificacion':verify}


def _interactivo(out,perfiles,kp,resumen,curves,usable,nombres,limite_tau,limite_altura,mostrar_todas):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from matplotlib.colors import to_hex
    def color(rid):return to_hex(_color(rid))
    common=dict(template='plotly_white',height=850,margin=dict(l=70,r=20,t=80,b=65),
                legend=dict(orientation='h',y=-.15),hovermode='closest')
    f=make_subplots(rows=1,cols=2,specs=[[{'type':'xy'},{'type':'table'}]],column_widths=[.64,.36],horizontal_spacing=.05)
    for rid,p in perfiles.items():
        f.add_trace(go.Scatter(x=p.chi_rel_m,y=p.elevacion_m,mode='lines',name=f'R{rid:02d} · {nombres[rid]}',
             legendgroup=str(rid),line=dict(color=color(rid),width=1.5),hovertemplate='χ relativo %{x:.1f} m<br>z %{y:.1f} m<extra>%{fullData.name}</extra>'),1,1)
        m=kp[kp.rid.eq(rid)&kp.valido]
        if len(m):f.add_trace(go.Scatter(x=m.chi_rel_KP_m,y=m.z_KP_m,mode='markers+text',text=m.rotulo,
             textposition='top center',customdata=m[['id_kp']].to_numpy(),legendgroup=str(rid),showlegend=False,
             marker=dict(symbol='star',size=14,color=color(rid),line=dict(color='black',width=1)),
             hovertemplate='%{customdata[0]}<br>χ relativo %{x:.1f} m<br>Altura %{y:.1f} m<extra></extra>'),1,1)
    rows=[]
    for _,v in resumen.iterrows():
        rows.append([v.rotulo,_nombre(v.rio),f'{v.z_KP_m:.1f}' if np.isfinite(v.z_KP_m) else '—',
          _rango([v.tau_5e_min_ka,v.tau_5e_max_ka]),_rango([v.tau_5c_min_ka,v.tau_5c_max_ka])])
    f.add_trace(go.Table(header=dict(values=['KP / río','Nombre','z (m)','τ 5e (ka)','τ 5c (ka)'],fill_color='#203d51',font_color='white'),
       cells=dict(values=list(map(list,zip(*rows))),height=31,fill_color='#f0f5f9',font_size=11),columnwidth=[55,125,55,90,90]),1,2)
    f.update_layout(**common,title='Perfiles χ · origen horizontal en cada outlet analítico; alturas sin desplazar')
    f.update_xaxes(title='χ − χ(outlet analítico) (m)',rangemode='tozero',row=1,col=1)
    f.update_yaxes(title='Elevación DEM (m)',range=[0,limite_altura] if limite_altura else None,row=1,col=1)
    t=go.Figure();seen=set()
    boundaries={}
    for rid,v,tau,z in curves:
        key=(rid,v.caso_B)
        boundaries.setdefault(key,[]).append((v.K_m02_anio,v.id_medicion))
    draw_ids={v[1] for values in boundaries.values() for v in [min(values),max(values)]}
    for rid,v,tau,z in curves:
        # Reducción exclusivamente gráfica; los KP se agregan por separado con coordenadas exactas.
        idx=np.unique(np.r_[np.arange(0,len(z),max(1,len(z)//1500)),len(z)-1])
        if not mostrar_todas and v.id_medicion not in draw_ids:continue
        key=str(rid);first=key not in seen;seen.add(key)
        t.add_trace(go.Scattergl(x=tau[idx],y=z[idx],mode='lines',name=f'R{rid:02d} · {nombres[rid]}',
            legendgroup=key,showlegend=first,opacity=.36,line=dict(color=color(rid),width=1,dash='solid' if v.mis=='5e' else 'dash'),
            hovertemplate=f'{v.id_medicion} · {v.id_B_corto} · MIS {v.mis}<br>K={v.K_m02_anio:.4e}<br>τ=%{{x:.1f}} ka<br>z=%{{y:.1f}} m<extra></extra>'))
    for (rid,mis),q in usable.groupby(['rid','mis'],sort=True):
        t.add_trace(go.Scatter(x=q.tau_KP_ka,y=q.z_kp_m,mode='markers',name=f'KP R{int(rid):02d} · MIS {mis}',
            legendgroup=str(int(rid)),showlegend=False,
            customdata=q[['id_kp','id_medicion','id_B','K_m02_anio']].to_numpy(),
            marker=dict(color=color(rid),size=12,symbol='star' if mis=='5e' else 'diamond',line=dict(color='black',width=.7)),
            hovertemplate='%{customdata[0]}<br>%{customdata[1]} · B %{customdata[2]}<br>K=%{customdata[3]:.4e}<br>τ=%{x:.1f} ka<br>z=%{y:.1f} m<extra></extra>'))
        lo,hi=q.tau_KP_ka.min(),q.tau_KP_ka.max()
        t.add_annotation(x=(lo+hi)/2,y=q.z_kp_m.iloc[0],text=f'R{int(rid):02d} · {mis}',showarrow=True,ay=-35,ax=20,bgcolor='white',bordercolor=color(rid),font_size=12)
    t.update_layout(**common,title='Perfiles τ · τ=0 en cada B; todos los escenarios K; doble clic para recuperar la vista')
    t.update_xaxes(title='τ relativo a cada B (ka)',range=[0,limite_tau] if limite_tau else None,rangemode='tozero')
    t.update_yaxes(title='Elevación DEM (m)',range=[0,limite_altura] if limite_altura else None,rangemode='tozero')
    config=dict(scrollZoom=True,displaylogo=False,toImageButtonOptions=dict(format='png',scale=2))
    html='''<!doctype html><html lang="es"><meta charset="utf-8"><title>KP · origen común</title>
<style>body{font:16px system-ui;margin:25px;color:#203d51}button{padding:12px;margin:6px;border:1px solid #54758d;background:#edf4fa;cursor:pointer}.panel{display:none}.panel.activo{display:block}p{max-width:1400px;line-height:1.6}</style>
<h1>Perfiles superpuestos y tiempos de propagación en los KP</h1>
<p>La elevación conserva el DEM original. χ se alinea en el outlet analítico de cada río; τ se alinea en cada B reconstruido. Estos orígenes coinciden numéricamente, no geográficamente. Los rangos son sensibilidad a los escenarios B/K; no intervalos de confianza ni edades tectónicas independientes.</p>
<button onclick="ver('chi')">18 ríos · alturas de los KP y tabla de τ</button><button onclick="ver('tau')">Ríos con K · perfiles τ y escenarios</button>
<p>Rueda: zoom · arrastrar: desplazar/seleccionar · clic en leyenda: mostrar/ocultar río · doble clic: aislar río o recuperar vista. En τ, pasar sobre un KP muestra el ID de medición y B.</p>'''
    html+='<div id="chi" class="panel activo">'+f.to_html(full_html=False,include_plotlyjs=True,config=config)+'</div>'
    html+='<div id="tau" class="panel">'+t.to_html(full_html=False,include_plotlyjs=False,config=config)+'</div>'
    html+='''<p>τ(KP) = [χ(KP) − χ(B)] / (K·A₀^0,40) / 1000, en ka; n=1 y A₀=10⁶ m². El MIS de la tabla identifica la terraza usada para calibrar K, no una asignación del KP a ese MIS. «—» significa que no hay τ disponible. Chanco no tiene KP manual.</p>
<script>function ver(id){document.querySelectorAll('.panel').forEach(e=>e.classList.toggle('activo',e.id===id));document.querySelectorAll('#'+id+' .plotly-graph-div').forEach(e=>Plotly.Plots.resize(e));}</script></html>'''
    target=out/'Perfiles_KP_origen_comun_interactivo.html';target.write_text(html,encoding='utf-8')
    return target


def enlazar_panel(salida):
    """Añade el gráfico al índice existente sin sustituir resultados anteriores."""
    salida=Path(salida);page=salida/'ABRIR_RESULTADOS.html'
    out=salida/'Comparacion_KP_origen_comun'
    if not page.exists() or not (out/'Perfiles_chi_origen_comun_KP.png').exists():return
    text=page.read_text(encoding='utf-8')
    if '<!-- KP_ORIGEN_COMUN -->' in text:return
    section='''<!-- KP_ORIGEN_COMUN --><h2>Perfiles con origen común y KP rotulados</h2>
<p><a href="Comparacion_KP_origen_comun/Perfiles_KP_origen_comun_interactivo.html">Abrir comparación interactiva de perfiles y tiempos en los KP</a></p>
<img src="Comparacion_KP_origen_comun/Perfiles_chi_origen_comun_KP.png">
<img src="Comparacion_KP_origen_comun/Perfiles_tau_origen_B_KP.png">'''
    pos=text.find('<img ')
    text=text[:pos]+section+'\n'+text[pos:] if pos>=0 else text+'\n'+section
    page.write_text(text,encoding='utf-8')
