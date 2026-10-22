"""Visualizar datos ya calculados, sin ajustar pendientes ni recalcular chi.

Importar desde Perfiles_y_fallas_20261013.ipynb. No escribe archivos al importar.
"""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import base64
import io
import numpy as np
import pandas as pd


def preparar(sw, perfiles, segmentos, theta, ventana_m, percentiles=(2, 100), limites=None):
    sw.validar_ksn(perfiles, ventana_m, theta)
    pct = np.asarray(percentiles, dtype=float)
    if len(pct) != 2 or not 0 <= pct[0] < pct[1] <= 100:
        raise ValueError('Usar dos percentiles crecientes entre 0 y 100.')
    valores = np.concatenate([d['ksn_040'] for d in segmentos.values()])
    valores = valores[np.isfinite(valores) & (valores > 0)]
    if not len(valores):
        raise ValueError('No hay ksn positivos disponibles.')
    vmin, vmax = np.percentile(valores, pct) if limites is None else np.asarray(limites, float)
    if not np.isfinite([vmin, vmax]).all() or vmin <= 0 or vmax < vmin:
        raise ValueError('Los límites de color deben ser positivos, finitos y crecientes.')
    if vmax == vmin:
        vmin, vmax = vmin / 2, vmax * 2
    return dict(theta=float(theta), ventana_m=float(ventana_m), vmin=float(vmin),
                vmax=float(vmax), minimo_real=float(valores.min()), maximo_real=float(valores.max()))


def _ids(ids, datos):
    seleccion = sorted(datos, key=int) if ids is None else list(map(str, ids))
    if not seleccion or set(seleccion) - set(datos) or len(seleccion) != len(set(seleccion)):
        raise ValueError('Revisar selección de ríos: vacía, repetida o no disponible.')
    return seleccion


def _validar(vista, perfiles, sw):
    sw.validar_ksn(perfiles, vista['ventana_m'], vista['theta'])


def _partes(geom):
    if geom.is_empty:
        return
    if geom.geom_type == 'LineString':
        yield geom
    elif hasattr(geom, 'geoms'):
        for part in geom.geoms:
            yield from _partes(part)


def _mascaras_segmentos(d, kps):
    mitad = (np.asarray(d['s0_m']) + np.asarray(d['s1_m'])) / 2
    if kps.empty:
        return np.zeros(len(mitad), bool)
    # Puede haber varias selecciones; la máscara solo sirve para indicar contexto.
    return mitad < float(kps.distancia_m.max())


def _cruces(perfil, trazas, maule):
    """Cruces geométricos XY; no asignan falla, bloque ni origen tectónico."""
    from shapely.geometry import LineString
    xy = perfil[['x_utm', 'y_utm']].to_numpy(float)
    linea = LineString(xy)
    cumul = np.r_[0., np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))]
    s = perfil.distancia_m.to_numpy(float)
    registros = []
    def puntos(geom):
        if geom.is_empty:
            return
        if geom.geom_type == 'Point':
            yield geom
        elif hasattr(geom, 'geoms'):
            for part in geom.geoms:
                yield from puntos(part)
        # Una superposición lineal no tiene un punto de cruce único.
    for capa, contexto in [(trazas, 'CHAF'), (maule, 'Maule')]:
        if capa is None:
            continue
        for r in capa.itertuples():
            etiqueta = str(getattr(r, 'familia', getattr(r, 'falla', contexto)))
            for pt in puntos(linea.intersection(r.geometry)):
                distancia = float(np.interp(linea.project(pt), cumul, s))
                if not any(abs(distancia - d) < .1 and etiqueta == e for d, e in registros):
                    registros.append((distancia, etiqueta))
    return sorted(registros)


def _marcar_cruces(axes, cruces, perfil=None, chi=None):
    for i, (s, etiqueta) in enumerate(cruces):
        pos = s / 1000 if chi is None else np.interp(s, perfil.distancia_m, chi)
        for ax in axes:
            ax.axvline(pos, color='#783ba1', ls=':', lw=.9, alpha=.75)
        axes[0].text(pos, .98 - (i % 3) * .09, etiqueta,
                     transform=axes[0].get_xaxis_transform(), va='top', ha='left',
                     fontsize=7, rotation=90, color='#783ba1', clip_on=True)


def mapa(sw, vista, perfiles, segmentos, manuales, nombres, trazas, maule,
         ruta_dem, ids=None, solo_bajo_kp=False, ancho=5):
    import plotly.graph_objects as go
    from plotly.colors import sample_colorscale
    _validar(vista, perfiles, sw)
    seleccion = _ids(ids, perfiles)
    fig = sw.mapa_ksn_interactivo(perfiles_ksn=perfiles, segmentos_ksn=segmentos,
        manuales=manuales, nombres=nombres, fallas_chaf=trazas, fallas_maule=maule,
        ruta_dem=ruta_dem, ventana_m=vista['ventana_m'], theta=vista['theta'], rios=seleccion)
    # Conservar hillshade, cauces, KP y geometrías. Sustituir solo las trazas de color.
    fig.data = tuple(t for t in fig.data if t.customdata is None and not t.marker.showscale
                     and t.name != 'KP manuales')
    if solo_bajo_kp:
        fig.data[0].name = 'Cauce completo; gris = fuera de selección o sin ksn válido'
    vistos = set()
    for t in fig.data:
        nombre = t.name or ''
        if '| CHAF:' in nombre:
            id_traza = nombre.split('|', 1)[1].strip()
            registro = trazas.loc[trazas.id_traza.eq(id_traza)]
            familia = nombre.split('|')[0].strip().replace('*', '')
            buz = str(registro.iloc[0].get('dir_buz', '')).strip().upper() if len(registro) else ''
            grupo = f'{familia} ({buz})' if buz else familia
            t.hovertemplate = nombre + '<extra></extra>'
            t.name, t.legendgroup = grupo, grupo
            t.showlegend = grupo not in vistos
            vistos.add(grupo)
            t.line.color = '#b83b35' if familia == 'Carranza' and buz != 'NW' else '#a17b00' if familia == 'Carranza' else '#783ba1'
        elif nombre == 'Maule: contexto':
            t.legendgroup = nombre
            t.showlegend = nombre not in vistos
            vistos.add(nombre)
    ncol = 48
    colores = sample_colorscale('Turbo', (np.arange(ncol) + .5) / ncol)
    v0, v1 = np.log10([vista['vmin'], vista['vmax']])
    bloques, valores, metas = [], [], []
    for rid in seleccion:
        d = segmentos[rid]
        k = np.asarray(d['ksn_040'], float)
        kps = manuales.loc[manuales.ID.eq(rid)]
        bajo = _mascaras_segmentos(d, kps)
        ok = np.isfinite(k) & (k > 0)
        if solo_bajo_kp:
            ok &= bajo
        if not ok.any():
            continue
        bloques.append(d['xy'][ok]); valores.append(k[ok])
        metas.append(np.column_stack([
            np.full(ok.sum(), nombres[rid], object), k[ok], d['s0_m'][ok], d['s1_m'][ok],
            np.where(d['ventana_parcial'][ok], 'sí', 'no'),
            np.where(d.get('ventana_cruza_KP', np.zeros(len(k), bool))[ok], 'sí', 'no')]))
    if not bloques:
        raise ValueError('No hay segmentos válidos para la selección; revisar KP y ríos.')
    xy, k, meta = np.concatenate(bloques), np.concatenate(valores), np.concatenate(metas)
    bins = np.clip(np.floor((np.log10(k) - v0) / (v1 - v0) * ncol).astype(int), 0, ncol - 1)
    for b, color in enumerate(colores):
        sel = bins == b
        if not sel.any():
            continue
        coords = np.concatenate([xy[sel], np.full((sel.sum(), 1, 2), np.nan)], axis=1).reshape(-1, 2) / 1000
        fig.add_trace(go.Scattergl(x=coords[:, 0], y=coords[:, 1], mode='lines',
            line=dict(color=color, width=ancho), customdata=np.repeat(meta[sel], 3, axis=0),
            showlegend=False, connectgaps=False,
            hovertemplate='%{customdata[0]}<br>ksn: %{customdata[1]:.3g}<br>'
                'Distancia desde outlet: %{customdata[2]:.0f}–%{customdata[3]:.0f} m<br>'
                'Ventana parcial: %{customdata[4]}<br>Ventana cruza KP: %{customdata[5]}<extra></extra>'))
    # KP por encima de los cauces coloreados.
    kp = manuales.loc[manuales.ID.isin(seleccion)]
    if not kp.empty:
        fig.add_trace(go.Scattergl(x=kp.x_utm / 1000, y=kp.y_utm / 1000, mode='markers+text',
            text=kp.ID, textposition='top right', marker=dict(symbol='star', size=12, color='white',
            line=dict(color='black', width=1)), name='KP manuales', hovertext=kp.id_kp, hoverinfo='text'))
    ticks = np.geomspace(vista['vmin'], vista['vmax'], 6)
    labels = [f'{v:.3g}' for v in ticks]
    if vista['vmin'] > vista['minimo_real']:
        labels[0] = '≤' + labels[0]
    if vista['vmax'] < vista['maximo_real']:
        labels[-1] = '≥' + labels[-1]
    punto = xy[0, 0] / 1000
    fig.add_trace(go.Scattergl(x=[punto[0]] * 2, y=[punto[1]] * 2, mode='markers',
        marker=dict(size=1, opacity=0, color=[v0, v1], colorscale='Turbo', cmin=v0, cmax=v1,
            showscale=True, colorbar=dict(title=f"ksn (m^{2*vista['theta']:g})<br>escala logarítmica",
                tickvals=np.log10(ticks), ticktext=labels, x=1.02, y=.5, len=.72)),
        showlegend=False, hoverinfo='skip'))
    texto = 'color solo aguas abajo de KP' if solo_bajo_kp else 'cauce completo'
    fig.update_layout(title=f"ksn local | {texto} | ventana {vista['ventana_m']:g} m | theta {vista['theta']:g}",
        legend=dict(orientation='h', x=0, y=-.13, yanchor='top', font=dict(size=10)),
        margin=dict(l=70, r=140, t=100, b=150), height=1050)
    return fig


def _fondo(fig_mapa):
    from PIL import Image
    if not fig_mapa.layout.images:
        return None
    im = fig_mapa.layout.images[0]
    if not isinstance(im.source, str) or not im.source.startswith('data:image'):
        return None
    img = np.asarray(Image.open(io.BytesIO(base64.b64decode(im.source.split(',', 1)[1]))))
    return img, (float(im.x), float(im.x + im.sizex), float(im.y - im.sizey), float(im.y)), float(im.opacity)


def _figura_guardar(vista, fig, nombre):
    vista.setdefault('figuras', []).append((nombre, fig))


def paneles_planta(vista, fig_mapa, perfiles, segmentos, manuales, nombres, trazas, maule,
                   ids=None, solo_bajo_kp=False):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable
    seleccion = _ids(ids, perfiles)
    fondo = _fondo(fig_mapa)
    norm = LogNorm(vista['vmin'], vista['vmax'])
    for pagina, inicio in enumerate(range(0, len(seleccion), 6), 1):
        lote = seleccion[inicio:inicio + 6]
        fig, axes = plt.subplots(int(np.ceil(len(lote) / 2)), 2, figsize=(13, 4.8 * int(np.ceil(len(lote) / 2))),
                                 squeeze=False, constrained_layout=True)
        for ax, rid in zip(axes.ravel(), lote):
            p, d = perfiles[rid], segmentos[rid]
            xy = p[['x_utm', 'y_utm']].to_numpy(float) / 1000
            if fondo is not None:
                ax.imshow(fondo[0], extent=fondo[1], origin='upper', alpha=fondo[2], zorder=0)
            ax.plot(xy[:, 0], xy[:, 1], color='.65', lw=1, zorder=1)
            k = np.asarray(d['ksn_040'])
            kp = manuales.loc[manuales.ID.eq(rid)]
            ok = np.isfinite(k) & (k > 0)
            if solo_bajo_kp:
                ok &= _mascaras_segmentos(d, kp)
            coleccion = LineCollection(d['xy'][ok] / 1000, cmap='turbo', norm=norm, linewidths=2.5, zorder=3)
            coleccion.set_array(k[ok]); ax.add_collection(coleccion)
            for capa, color in [(trazas, '#783ba1'), (maule, '.4')]:
                if capa is None:
                    continue
                for geom in capa.geometry:
                    for part in _partes(geom):
                        v = np.asarray(part.coords)[:, :2] / 1000
                        ax.plot(v[:, 0], v[:, 1], '--', color=color, lw=1, zorder=2)
            if not kp.empty:
                ax.scatter(kp.x_utm / 1000, kp.y_utm / 1000, marker='*', s=110, c='white', edgecolors='black', zorder=4)
            pad = max(.5, .06 * np.ptp(xy, axis=0).max())
            ax.set(xlim=(xy[:, 0].min() - pad, xy[:, 0].max() + pad),
                   ylim=(xy[:, 1].min() - pad, xy[:, 1].max() + pad),
                   title=f'{rid}: {nombres[rid]}' + (' | sin KP' if kp.empty else ''),
                   xlabel='Este UTM (km)', ylabel='Norte UTM (km)')
            ax.set_aspect('equal'); ax.grid(alpha=.15); ax.ticklabel_format(useOffset=False, style='plain')
        for ax in axes.ravel()[len(lote):]:
            ax.set_visible(False)
        cb = fig.colorbar(ScalarMappable(norm=norm, cmap='turbo'), ax=list(axes.ravel()[:len(lote)]), shrink=.8)
        cb.set_label(f"ksn local (m^{2*vista['theta']:g}) — escala común")
        fig.suptitle(f"Mapa por río | ventana {vista['ventana_m']:g} m | página {pagina}")
        _figura_guardar(vista, fig, f'mapas_por_rio_{pagina:02d}')
        plt.show(); plt.close(fig)


def _marcar_kp(axes, kp, largo_m, ventana_m):
    for r in kp.itertuples():
        for ax in axes:
            ax.axvline(r.distancia_m / 1000, color='black', ls='--', lw=.9)
            ax.axvspan(max(0, r.distancia_m - ventana_m / 2) / 1000,
                       min(largo_m, r.distancia_m + ventana_m / 2) / 1000, color='gray', alpha=.12)


def paneles_perfiles(vista, perfiles, manuales, nombres, ids=None, trazas=None, maule=None):
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LogNorm
    seleccion = _ids(ids, perfiles)
    norm = LogNorm(vista['vmin'], vista['vmax'])
    for pagina, inicio in enumerate(range(0, len(seleccion), 4), 1):
        lote = seleccion[inicio:inicio + 4]
        fig, axes = plt.subplots(len(lote), 2, figsize=(14, 3.4 * len(lote)), squeeze=False, constrained_layout=True)
        for fila, rid in enumerate(lote):
            q = perfiles[rid]
            x, z, k = q.distancia_m.to_numpy(float) / 1000, q.elevacion_m.to_numpy(float), q.ksn_040.to_numpy(float)
            a, b = axes[fila]
            a.plot(x, z, color='.7', lw=1)
            coords = np.column_stack([x, z]); pares = np.stack([coords[:-1], coords[1:]], axis=1)
            valores = (k[:-1] + k[1:]) / 2
            ok = np.isfinite(valores) & (valores > 0)
            lc = LineCollection(pares[ok], cmap='turbo', norm=norm, linewidths=2)
            lc.set_array(valores[ok]); a.add_collection(lc)
            b.plot(x, k, color='darkorange', lw=1.2)
            kp = manuales.loc[manuales.ID.eq(rid)]
            _marcar_kp([a, b], kp, q.distancia_m.iloc[-1], vista['ventana_m'])
            _marcar_cruces([a, b], _cruces(q, trazas, maule))
            a.set(title=f'{rid}: {nombres[rid]}', ylabel='Elevación (m)')
            b.set(title='ksn local sin recorte de valores', ylabel=f"ksn (m^{2*vista['theta']:g})")
            for ax in (a, b):
                ax.set(xlabel='Distancia desde outlet (km)', xlim=(x[0], x[-1])); ax.grid(alpha=.2)
        fig.suptitle('Negro: KP | gris: ventanas que cruzan KP | violeta: cruce XY con traza cartográfica')
        _figura_guardar(vista, fig, f'perfiles_completos_{pagina:02d}')
        plt.show(); plt.close(fig)


def detalle(vista, perfiles, manuales, nombres, ajustes, tabla, ids=None, margen_kp_m=2500,
            trazas=None, maule=None):
    import matplotlib.pyplot as plt
    seleccion = _ids(ids, perfiles)
    for rid in seleccion:
        p = perfiles[rid]; kp_rio = manuales.loc[manuales.ID.eq(rid)]
        if kp_rio.empty:
            print(f'{rid}: {nombres[rid]} sin KP; consultar su perfil completo.'); continue
        for kp in kp_rio.itertuples():
            clave = (rid, kp.id_kp)
            if clave not in ajustes:
                print(f'{clave}: sin ajuste; revisar estado en tabla.'); continue
            a = ajustes[clave]
            if not np.isclose(tabla.loc[tabla.ID.eq(rid) & tabla.id_kp.eq(kp.id_kp), 'theta'].iloc[0], vista['theta']):
                raise ValueError('Los ajustes tienen otro theta: recalcular sección 5.')
            s = p.distancia_m.to_numpy(float)
            mask = a['mask']; chi = a['chi']; z = p.elevacion_m.to_numpy(float)
            fig, axes = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)
            a0, a1, a2, a3 = axes.ravel()
            a0.plot(s / 1000, z, color='.65', label='Cauce')
            a0.plot(s[mask] / 1000, z[mask], color='steelblue', lw=2, label='Tramo ajustado')
            a0.scatter(kp.distancia_m / 1000, kp.elevacion_m, marker='*', s=110, color='crimson', label='KP')
            izquierda = max(s[0], kp.distancia_m - margen_kp_m)
            derecha = min(s[-1], kp.distancia_m + margen_kp_m)
            _marcar_kp([a0, a1], kp_rio, s[-1], vista['ventana_m'])
            cruces = _cruces(p, trazas, maule)
            _marcar_cruces([a0, a1], cruces)
            a0.set(xlim=(izquierda / 1000, derecha / 1000), ylabel='Elevación (m)', xlabel='Distancia desde outlet (km)', title='Zoom alrededor del KP')
            a0.legend(fontsize=8)
            a1.plot(s / 1000, p.ksn_040, color='darkorange', label='ksn local')
            t = tabla.loc[tabla.ID.eq(rid) & tabla.id_kp.eq(kp.id_kp)].iloc[0]
            if np.isfinite(t.get('ksn_040', np.nan)):
                a1.plot(s[mask] / 1000, np.full(mask.sum(), t.ksn_040), 'k:', label='ksn del ajuste χ–z')
            a1.set(xlim=(izquierda / 1000, derecha / 1000), ylabel=f"ksn (m^{2*vista['theta']:g})", xlabel='Distancia desde outlet (km)', title='Local y representativo son estimadores distintos')
            a1.legend(fontsize=8)
            a2.plot(chi[mask], z[mask], color='steelblue', label='Tramo')
            a2.plot(chi[mask], a['fit'].intercept + a['fit'].slope * chi[mask], 'k--', label='Ajuste')
            a2.scatter(chi[a['j']], z[a['j']], marker='*', s=110, color='crimson')
            a2.set(xlabel='χ (m)', ylabel='Elevación (m)', title=f"ksn tramo={t.ksn_040:.2f} | R²={t.R2:.3f} | RMSE={t.RMSE_m:.2f} m")
            a2.legend(fontsize=8)
            a3.plot(chi[mask], a['residuos'], color='steelblue')
            a3.axhline(0, color='black', ls='--', lw=.8)
            a3.set(xlabel='χ (m)', ylabel='Residuo (m)', title='Revisar estructura de residuos, no solo R²')
            cruces_tramo = [(d, e) for d, e in cruces if s[mask].min() <= d <= s[mask].max()]
            _marcar_cruces([a2, a3], cruces_tramo, perfil=p, chi=chi)
            for ax in axes.ravel():
                ax.grid(alpha=.2)
            fig.suptitle(f'{rid}: {nombres[rid]} | {kp.id_kp} | ventana {vista["ventana_m"]:g} m')
            _figura_guardar(vista, fig, f'detalle_R{int(rid):02d}_{kp.id_kp}')
            plt.show(); plt.close(fig)


def comparar(vista, tabla, ids=None):
    import matplotlib.pyplot as plt
    if 'ksn_040' not in tabla:
        print('No hay ajustes disponibles.'); return
    q = tabla.loc[np.isfinite(tabla.ksn_040) & tabla.ksn_040.gt(0)].copy()
    if ids is not None:
        q = q.loc[q.ID.isin(list(map(str, ids)))]
    if q.empty:
        print('No hay ajustes positivos para la selección.'); return
    if not np.allclose(q.theta, vista['theta']):
        raise ValueError('Tabla y mapa utilizan theta diferente.')
    q['orden'] = pd.to_numeric(q.ID)
    q = q.sort_values(['orden', 'id_kp'])
    y = np.arange(len(q))
    etiquetas = q.ID + ': ' + q.rio
    repetidos = q.ID.duplicated(keep=False)
    etiquetas = etiquetas.where(~repetidos, etiquetas + ' | ' + q.id_kp)
    fig, axes = plt.subplots(1, 2, figsize=(14, max(5, .42 * len(q))), sharey=True, constrained_layout=True)
    axes[0].barh(y, q.ksn_040, color='steelblue')
    axes[1].barh(y, q.longitud_ajuste_m / 1000, color='gray')
    axes[0].set_yticks(y, etiquetas); axes[0].invert_yaxis()
    axes[0].set(xlabel=f"ksn del tramo (m^{2*vista['theta']:g})", title='Resumen por tramo bajo KP')
    axes[1].set(xlabel='Longitud ajustada (km)', title='Extensión que representa cada valor')
    for i, r in enumerate(q.itertuples()):
        axes[0].text(r.ksn_040, i, f' {r.ksn_040:.2f}', va='center', fontsize=9)
    axes[0].set_xlim(0, q.ksn_040.max() * 1.16)
    for ax in axes:
        ax.grid(axis='x', alpha=.2)
    fig.suptitle('Valores exploratorios; las longitudes y los perfiles deben revisarse antes de interpretar')
    _figura_guardar(vista, fig, 'comparacion_tramos')
    plt.show(); plt.close(fig)


def rumbo(vista, comparacion, alcances):
    import matplotlib.pyplot as plt
    if comparacion.empty:
        print('Figura por rumbo pendiente: completar asociaciones revisadas y ejecutar su cálculo.'); return
    familias = list(comparacion.falla.unique())
    fig, axes = plt.subplots(1, len(familias), figsize=(7 * len(familias), 5), squeeze=False,
                             sharey=True, constrained_layout=True)
    for ax, falla in zip(axes.ravel(), familias):
        q = comparacion.loc[comparacion.falla.eq(falla)]
        if not np.allclose(q.theta, vista['theta']):
            raise ValueError('Comparación por rumbo calculada con otro theta.')
        for familia, lo, hi in alcances:
            if familia == falla:
                ax.axvspan(lo, hi, color='gray', alpha=.08)
        ax.scatter(q.distancia_rumbo_km, q.ksn_040, color='steelblue')
        for r in q.itertuples():
            ax.annotate(f'{r.ID}: {r.rio}', (r.distancia_rumbo_km, r.ksn_040),
                        xytext=(4, 4), textcoords='offset points', fontsize=8)
        ax.set(title=falla, xlabel='Distancia proyectada sobre rumbo (km)', ylabel=f"ksn del tramo (m^{2*vista['theta']:g})")
        ax.grid(alpha=.2)
    fig.suptitle('Sombreado: registros cartográficos disponibles; no longitud total de falla')
    _figura_guardar(vista, fig, 'comparacion_rumbo')
    plt.show(); plt.close(fig)


def guardar(vista, carpeta, mapas):
    from datetime import datetime
    import json
    import memoria_comun as pa
    salida = pa.nueva_carpeta(carpeta)
    salida.mkdir(parents=True, exist_ok=True)
    for nombre, fig in mapas.items():
        fig.write_html(salida / f'{nombre}_{salida.name}.html', include_plotlyjs=True)
    for nombre, fig in vista.get('figuras', []):
        fig.savefig(salida / f'{nombre}_{salida.name}.png', dpi=300, bbox_inches='tight')
        fig.savefig(salida / f'{nombre}_{salida.name}.pdf', bbox_inches='tight')
    parametros = {k: v for k, v in vista.items() if k != 'figuras'}
    (salida / 'parametros_visualizacion.json').write_text(json.dumps(parametros, indent=2), encoding='utf-8')
    print('Figuras guardadas en:', salida)
    return salida
