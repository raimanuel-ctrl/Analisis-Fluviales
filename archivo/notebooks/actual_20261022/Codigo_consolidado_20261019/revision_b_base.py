"""Revisión geomorfológica de B y escenarios K; no selecciona paleocostas automáticamente."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.windows import from_bounds
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import map_coordinates
from scipy.stats import linregress
from shapely.geometry import LineString, Point
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib.ticker import ScalarFormatter

CRS = 'EPSG:32718'
THETA = .40
A0 = 1_000_000.


def huella(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def mis_normal(value):
    m = re.search(r'5\s*([ec])', str(value).lower())
    if not m:
        raise ValueError(f'MIS no reconocido: {value!r}; usar 5e o 5c.')
    return '5' + m.group(1)


def bandera(value):
    return str(value).strip().lower() in {'1', '1.0', 'true', 'si', 'sí', 'yes'}


def cargar_perfiles(base, ids):
    base = Path(base)
    cfg = json.loads((base / 'Configuracion.json').read_text(encoding='utf-8-sig'))
    mfpath = base / 'Productos' / 'Datos' / cfg['run_id'] / 'perfiles' / f"manifest_{cfg['run_id']}.json"
    mf = json.loads(mfpath.read_text(encoding='utf-8-sig'))
    if mf.get('estado') != 'completo':
        raise ValueError('El manifiesto de perfiles no está completo.')
    for k, target in [('theta', THETA), ('A0_m2', A0)]:
        if not np.isclose(mf['config'][k], target):
            raise ValueError(f'{k} no coincide con los perfiles guardados.')
    profiles, names, paths = {}, {}, {}
    for rid in ids:
        item = mf['files'][f'rio_{int(rid)}']
        path = Path(item['path'])
        if not path.is_absolute():
            path = mfpath.parent / path
        if huella(path) != item['sha256']:
            raise ValueError(f'Cambió el perfil R{rid:02d}: {path}')
        p = preparar_perfil(pd.read_csv(path, float_precision='round_trip'))
        profiles[int(rid)] = p
        names[int(rid)] = mf['rios'][str(int(rid))]['nombre']
        paths[int(rid)] = path
    return profiles, names, paths, cfg, mfpath


def preparar_perfil(p):
    p = p.copy()
    cols = ['distancia_m', 'elevacion_m', 'area_m2', 'x_utm', 'y_utm']
    if len(p) < 2 or not np.isfinite(p[cols].to_numpy(float)).all():
        raise ValueError('Perfil incompleto: no se conectan nodos a través de huecos.')
    if not (np.diff(p.distancia_m) > 0).all() or (p.area_m2 <= 0).any():
        raise ValueError('Distancias no crecientes o área no positiva.')
    p['chi_m'] = A0 ** THETA * cumulative_trapezoid(p.area_m2.to_numpy() ** -THETA,
                                                   p.distancia_m, initial=0)
    p['s_geom_m'] = np.r_[0., np.cumsum(np.hypot(np.diff(p.x_utm), np.diff(p.y_utm)))]
    if not (np.diff(p.s_geom_m) > 0).all():
        raise ValueError('Nodos XY repetidos consecutivos.')
    return p


def linea_rio(p):
    return LineString(p[['x_utm', 'y_utm']].to_numpy())


def puntos_interseccion(g):
    if g.is_empty:
        return []
    if g.geom_type == 'Point':
        return [g]
    if g.geom_type in {'MultiPoint', 'GeometryCollection'}:
        out = []
        for sub in g.geoms:
            out.extend(puntos_interseccion(sub))
        return out
    raise ValueError('La marca se superpone con el cauce: se necesita un cruce transversal inequívoco.')


def capas_vacias():
    schemas = {
        'bordes_observados': {'rid': 'int', 'mis': 'str', 'grupo': 'str', 'margen': 'str',
                              'fuente': 'str', 'evidencia': 'str', 'revisado': 'int'},
        'prolongaciones_B': {'rid': 'int', 'mis': 'str', 'grupo': 'str', 'escenario': 'str',
                             'fuente': 'str', 'evidencia': 'str', 'revisado': 'int'},
        'transectos_revision': {'rid': 'int', 'mis': 'str', 'grupo': 'str', 'id_perfil': 'str',
                                'margen': 'str', 'ancho_m': 'float', 'plat_ini': 'float',
                                'plat_fin': 'float', 'acan_ini': 'float', 'acan_fin': 'float',
                                'nota': 'str'},
    }
    for layer, fields in schemas.items():
        dtype = {'int': 'int64', 'str': 'object', 'float': 'float64'}
        yield layer, gpd.GeoDataFrame({k: pd.Series(dtype=dtype[v]) for k, v in fields.items()},
                                     geometry=gpd.GeoSeries([], crs=CRS), crs=CRS)


def preparar_capas(folder, perfiles, nombres, ruta_marcas, terrazas, muestras, jm):
    """Idempotente: nunca sobreescribe una capa de digitización ya existente."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    editable = folder / 'Revision_B_EDITAR_EN_QGIS.gpkg'
    if not editable.exists():
        for layer, gdf in capas_vacias():
            gdf.to_file(editable, layer=layer, driver='GPKG', geometry_type='LineString')
    refpath = folder / 'Referencias_NO_EDITAR.gpkg'
    from pyogrio import list_layers
    existing = set(list_layers(refpath)[:, 0]) if refpath.exists() else set()
    original = gpd.read_file(ruta_marcas)
    if original.crs is None:
        raise ValueError('AnalisisK carece de CRS.')
    original = original.to_crs(CRS)
    original['fid_original'] = original.index.astype(int)
    original['tipo_B'] = 'cota; interpretacion inicial; no borde confirmado'
    # SQLite no distingue MIS/mis en nombres de campos: exportar una tabla normalizada.
    mref = muestras[['sample_id', 'mis', 'x_control', 'y_control', 'AGE (ka)', 'error', 'assigned MIS']].copy()
    references = {
        'rios': gpd.GeoDataFrame([{'rid': r, 'rio': nombres[r], 'geometry': linea_rio(p)}
                                 for r, p in perfiles.items()], crs=CRS),
        'B_inicial_por_cota': original,
        'Terrazas_JM': gpd.GeoDataFrame(terrazas, geometry=gpd.points_from_xy(
            terrazas.este_utm18s_m, terrazas.norte_utm18s_m), crs=CRS),
        'Muestras_compatibles': gpd.GeoDataFrame(mref, geometry=gpd.points_from_xy(
            mref.x_control, mref.y_control), crs=CRS),
        'JM_estacion_perfil_pagina': gpd.GeoDataFrame(jm, geometry=gpd.points_from_xy(jm.x, jm.y), crs=CRS)}
    for layer, frame in references.items():
        # También completa preparaciones interrumpidas sin tocar capas ya terminadas.
        if layer not in existing or len(gpd.read_file(refpath, layer=layer)) == 0:
            frame.to_file(refpath, layer=layer, driver='GPKG')
    return editable, refpath


def leer_capa(path, layer):
    g = gpd.read_file(path, layer=layer)
    if g.crs is None:
        raise ValueError(f'Falta CRS en {layer}.')
    g = g.to_crs(CRS)
    if g.empty:
        return g
    if g.geometry.isna().any() or g.geometry.is_empty.any() or not g.is_valid.all():
        raise ValueError(f'Geometrías vacías o inválidas en {layer}.')
    rid = pd.to_numeric(g.rid, errors='raise')
    if not np.isfinite(rid).all() or not np.equal(rid, np.round(rid)).all():
        raise ValueError(f'IDs no enteros en {layer}.')
    g['rid'] = rid.astype(int)
    g['mis'] = g.mis.map(mis_normal)
    for k in ['grupo', 'id_perfil', 'escenario']:
        if k in g and (g[k].isna().any() or g[k].astype(str).str.strip().eq('').any()):
            raise ValueError(f'Completar {k} en {layer}.')
    if not g.geom_type.eq('LineString').all():
        raise ValueError(f'{layer}: usar líneas simples, no multipartes.')
    return g


def leer_muestras(path):
    m = pd.read_csv(path, sep=';', decimal=',')
    m.columns = m.columns.str.strip()
    m['sample_id'] = m['Sample ID'].astype(str).str.strip()
    m['mis'] = m['assigned MIS'].map(mis_normal)
    # Verificación previa: los títulos Ny/Ex están invertidos respecto a sus valores.
    m['x_control'] = pd.to_numeric(m.Ny, errors='raise')
    m['y_control'] = pd.to_numeric(m.Ex, errors='raise')
    if not m.x_control.between(100000, 900000).all() or not m.y_control.between(5000000, 7000000).all():
        raise ValueError('Coordenadas de muestras incompatibles con UTM18S; revisar Ny/Ex.')
    if m.sample_id.duplicated().any():
        raise ValueError('Muestras duplicadas.')
    return m


def construir_B(perfiles, marcas, nuevas, incluir_inicial=True, excluir_fid=()):
    records, qc = [], []
    candidates = []
    if incluir_inicial:
        if marcas.crs is None:
            raise ValueError('Marcas iniciales sin CRS.')
        for fid, row in marcas.to_crs(CRS).iterrows():
            if int(fid) in excluir_fid:
                qc.append({'fid': int(fid), 'estado': 'excluido: otro experimento por decisión del usuario'})
                continue
            rid = int(row['id'])
            if float(row['id']) != rid:
                raise ValueError('ID original no entero.')
            candidates.append({'rid': rid, 'mis': mis_normal(row.MIS), 'grupo': 'inicial',
                'escenario': f'COTA_F{fid:02d}', 'fid_fuente': int(fid), 'geometry': row.geometry,
                'tipo_B': 'cota_inicial', 'revisado': False, 'fuente': 'AnalisisK.shp',
                'evidencia': 'prolongacion de cota entre terrazas; borde pendiente'})
    if nuevas is not None and not nuevas.empty:
        keys = ['rid', 'mis', 'grupo', 'escenario']
        if nuevas.duplicated(keys).any():
            raise ValueError('Prolongaciones repetidas: rid, MIS, grupo y escenario deben ser únicos.')
        for fid, row in nuevas.iterrows():
            candidates.append({**row.to_dict(), 'fid_fuente': int(fid), 'tipo_B': 'borde_reconstruido',
                               'revisado': bandera(row.revisado)})
    for row in candidates:
        rid = row['rid']
        if rid not in perfiles:
            qc.append({'rid': rid, 'estado': 'río no incluido en RIOS_PROCESAR'})
            continue
        river = linea_rio(perfiles[rid])
        hits = sorted(puntos_interseccion(river.intersection(row['geometry'])), key=river.project)
        if not hits:
            qc.append({'rid': rid, 'mis': row['mis'], 'escenario': row['escenario'],
                       'estado': 'sin cruce; no se ajusta ni se extrapola el cauce',
                       'distancia_al_cauce_m': row['geometry'].distance(river)})
        for j, hit in enumerate(hits, 1):
            p = perfiles[rid]
            along = river.project(hit)
            if row['revisado'] and not str(row.get('evidencia', '')).strip():
                raise ValueError('Un B marcado como revisado necesita una justificación en evidencia.')
            values = {key: float(np.interp(along, p.s_geom_m, p[col])) for key, col in
                      [('s_B_m', 'distancia_m'), ('chi_B_m', 'chi_m'), ('z_B_m', 'elevacion_m')]}
            label = f"R{rid:02d}_{row['mis']}_{row['grupo']}_{row['escenario']}_X{j:02d}"
            records.append({k: v for k, v in row.items() if k != 'geometry'} | values |
                           {'caso_B': label, 'x_B': hit.x, 'y_B': hit.y, 'n_cruces_linea': len(hits)})
    return pd.DataFrame(records), pd.DataFrame(qc)


def elegir_referencias(btab, muestras, terrazas, jm, manual_muestras=None, manual_terrazas=None):
    manual_muestras = manual_muestras or {}
    manual_terrazas = manual_terrazas or {}
    t = terrazas.copy()
    t['mis'] = t.mis_asignado.map(mis_normal)
    out = []
    for _, b in btab.iterrows():
        same_m = muestras[muestras.mis.eq(b.mis)].copy()
        same_t = t[t.mis.eq(b.mis)].copy()
        if same_m.empty or same_t.empty:
            raise ValueError(f'{b.caso_B}: faltan muestra o terraza del mismo MIS.')
        same_m['d'] = np.hypot(same_m.x_control-b.x_B, same_m.y_control-b.y_B)
        sm = manual_muestras.get(b.caso_B)
        choice = same_m[same_m.sample_id.eq(sm)] if sm else same_m.nsmallest(1, 'd')
        if len(choice) != 1:
            raise ValueError(f'Muestra manual inválida en {b.caso_B}.')
        m = choice.iloc[0]
        same_t['d'] = np.hypot(same_t.este_utm18s_m-m.x_control, same_t.norte_utm18s_m-m.y_control)
        st = manual_terrazas.get(b.caso_B)
        choice = same_t[same_t.terrace_id.eq(st)] if st else same_t.nsmallest(1, 'd')
        if len(choice) != 1:
            raise ValueError(f'Terraza manual inválida en {b.caso_B}.')
        v = choice.iloc[0]
        idx = jm[jm.terrace_id.eq(v.terrace_id)]
        source = idx.iloc[0] if len(idx) == 1 else pd.Series(dtype=object)
        vals = [v.elevacion_m, v.error_elevacion_m, v.edad_referencia_ka, v.incertidumbre_edad_ka]
        if not np.isfinite(vals).all() or vals[1] < 0 or vals[2] <= 0 or vals[3] < 0:
            raise ValueError('Altura, edad o magnitudes publicadas inválidas.')
        out.append(b.to_dict() | {'muestra_control': m.sample_id,
            'edad_sedimento_ka': m['AGE (ka)'], 'error_edad_sedimento_ka': m['error'],
            'x_muestra': m.x_control, 'y_muestra': m.y_control,
            'distancia_muestra_B_m': m['d'], 'terraza_id': v.terrace_id,
            'H_m': v.elevacion_m, 'error_H_publicado_m': v.error_elevacion_m,
            'T_ka': v.edad_referencia_ka, 'error_T_publicado_ka': v.incertidumbre_edad_ka,
            'x_terraza': v.este_utm18s_m, 'y_terraza': v.norte_utm18s_m,
            'distancia_terraza_muestra_m': v['d'],
            'seleccion_referencia': ('manual' if sm and st else 'proximidad con posibles reemplazos manuales'),
            'estacion_JM': source.get('station', ''), 'perfil_JM': source.get('profile', np.nan),
            'codigo_JM': source.get('shoreline_codigo', np.nan),
            'paginas_suplemento': source.get('paginas_pdf', ''),
            'fuente_H': v.fuente})
    return pd.DataFrame(out)


def cruces_cota(p, H):
    """Cruces exactos lineales. Las mesetas se conservan como intervalos ambiguos."""
    s, z, chi = [p[c].to_numpy(float) for c in ['distancia_m', 'elevacion_m', 'chi_m']]
    rows = []
    for j in range(len(p)-1):
        dz = z[j+1]-z[j]
        if dz == 0:
            if np.isclose(z[j], H, atol=1e-9, rtol=0):
                rows.append({'s_PT_m': s[j], 's_PT_fin_m': s[j+1], 'chi_PT_m': chi[j],
                    'chi_PT_fin_m': chi[j+1], 'direccion': 'meseta', 'Mchi_local': 0.,
                    'x_PT': p.x_utm.iloc[j], 'y_PT': p.y_utm.iloc[j], 'segmento': j})
            continue
        if min(z[j], z[j+1]) <= H <= max(z[j], z[j+1]):
            f = (H-z[j])/dz
            sp = s[j]+f*(s[j+1]-s[j])
            if rows and rows[-1]['direccion'] != 'meseta' and np.isclose(sp, rows[-1]['s_PT_m'], atol=1e-7, rtol=0):
                if dz > 0:
                    rows[-1]['direccion'] = 'ascendente'
                continue
            rows.append({'s_PT_m': sp, 's_PT_fin_m': sp, 'chi_PT_m': chi[j]+f*(chi[j+1]-chi[j]),
                'chi_PT_fin_m': chi[j]+f*(chi[j+1]-chi[j]),
                'x_PT': p.x_utm.iloc[j]+f*(p.x_utm.iloc[j+1]-p.x_utm.iloc[j]),
                'y_PT': p.y_utm.iloc[j]+f*(p.y_utm.iloc[j+1]-p.y_utm.iloc[j]),
                'direccion': 'ascendente' if dz > 0 else 'descendente',
                'Mchi_local': dz/(chi[j+1]-chi[j]), 'segmento': j})
    cols = ['cruce_PT', 's_PT_m', 's_PT_fin_m', 'chi_PT_m', 'chi_PT_fin_m',
            'x_PT', 'y_PT', 'direccion', 'Mchi_local', 'segmento']
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=cols)
    # Un nodo de meseta no debe pasar como punto inequívoco.
    flat = df[df.direccion.eq('meseta')]
    for _, a in flat.iterrows():
        touch = df.s_PT_m.between(a.s_PT_m-1e-7, a.s_PT_fin_m+1e-7) & df.direccion.ne('meseta')
        df = df[~touch]
    df = df.sort_values('s_PT_m').reset_index(drop=True)
    df.insert(0, 'cruce_PT', np.arange(1, len(df)+1))
    return df[cols]


def calcular_escenarios(refs, perfiles, seleccion_pt=None, datum_confirmado=False,
                        nota_datum='', nivel_mar=None):
    seleccion_pt = seleccion_pt or {}
    nivel_mar = nivel_mar or {'5e': {'valor_m': 5., 'error_m': 5.},
                             '5c': {'valor_m': -20., 'error_m': 5.}}
    rows, status = [], []
    if datum_confirmado and not str(nota_datum).strip():
        raise ValueError('Documentar evidencia de compatibilidad vertical antes de marcarla confirmada.')
    for _, b in refs.iterrows():
        cr = cruces_cota(perfiles[int(b.rid)], b.H_m)
        cr = cr[cr.s_PT_m > b.s_B_m+1e-6].copy()
        asc = cr[cr.direccion.eq('ascendente')]
        requested = seleccion_pt.get(b.caso_B)
        if requested is not None and requested not in asc.cruce_PT.tolist():
            raise ValueError(f'{b.caso_B}: P_T manual no es un cruce ascendente válido aguas arriba de B.')
        default = int(asc.cruce_PT.iloc[0]) if not asc.empty else None
        primary = requested if requested is not None else default
        status.append({'caso_B': b.caso_B, 'n_cruces_aguas_arriba': len(cr),
                       'n_ascendentes': len(asc), 'estado': 'con escenarios' if len(asc) else 'sin P_T ascendente puntual'})
        for _, pt in cr.iterrows():
            delta = pt.chi_PT_m-b.chi_B_m
            usable = pt.direccion != 'meseta' and delta > 0
            K = delta/(A0**THETA*b.T_ka*1000) if usable else np.nan
            approx_error = np.nan
            if usable and pt.Mchi_local != 0:
                echi = b.error_H_publicado_m/abs(pt.Mchi_local)
                approx_error = np.sqrt((echi/(A0**THETA*b.T_ka*1000))**2 +
                                       (K*b.error_T_publicado_ka/b.T_ka)**2)
            sea = nivel_mar[b.mis]
            U = (b.H_m-sea['valor_m'])/b.T_ka
            eU = np.sqrt((b.error_H_publicado_m/b.T_ka)**2 + (sea['error_m']/b.T_ka)**2 +
                        ((b.H_m-sea['valor_m'])*b.error_T_publicado_ka/b.T_ka**2)**2)
            q = perfiles[int(b.rid)]
            q = q[q.distancia_m.between(b.s_B_m, pt.s_PT_m)]
            ksn, r2, rmse = np.nan, np.nan, np.nan
            if len(q) >= 5 and q.chi_m.nunique() > 1:
                fit = linregress(q.chi_m, q.elevacion_m)
                if fit.slope > 0:
                    ksn, r2 = fit.slope*A0**THETA, fit.rvalue**2
                    rmse = np.sqrt(np.mean((q.elevacion_m-fit.intercept-fit.slope*q.chi_m)**2))
            principal = bool(pt.cruce_PT == primary and pt.direccion == 'ascendente')
            rows.append(b.to_dict() | pt.to_dict() | {'caso': f'{b.caso_B}_PT{int(pt.cruce_PT):03d}',
                'delta_chi_m': delta, 'K_George_m02_anio': K,
                'amplitud_K_propagada_aprox': approx_error, 'principal_operativo': principal,
                'seleccion_PT': 'manual' if requested is not None else 'primer ascendente provisional',
                'PT_elegible': pt.direccion == 'ascendente', 'datum_confirmado': bool(datum_confirmado),
                'nota_datum': nota_datum, 'estado_K': 'condicional; pendiente revision del modelo por propagacion',
                'U_terraza_mm_anio': U, 'amplitud_U_propagada_aprox': eU,
                'ksn_B_PT_m08': ksn, 'R2_B_PT': r2, 'RMSE_B_PT_m': rmse, 'n_fit': len(q),
                'U_Kksn_experimental_mm_anio': 1000*K*ksn if usable else np.nan,
                'nivel_mar_m': sea['valor_m'], 'amplitud_nivel_mar_m': sea['error_m']})
    return pd.DataFrame(rows), pd.DataFrame(status)


def sensibilidad_PT(res):
    rows = []
    if res.empty:
        return pd.DataFrame()
    for case, q in res[res.PT_elegible & res.K_George_m02_anio.gt(0)].groupby('caso_B'):
        k = q.K_George_m02_anio
        med = k.median()
        rows.append({'caso_B': case, 'n_PT': len(q), 'K_min': k.min(), 'K_mediana': med,
            'K_max': k.max(), 'rango_relativo_porcentaje': 100*(k.max()-k.min())/med,
            'separacion_PT_m': q.s_PT_m.max()-q.s_PT_m.min(),
            'criterio': 'rango de escenarios, no intervalo de confianza; sin aceptación automática'})
    return pd.DataFrame(rows)


def sensibilidad_B(refs, perfiles):
    """Aísla B: fija H,T y el mismo P_T para todos los B del río/MIS/referencia."""
    rows = []
    if refs.empty:
        return pd.DataFrame()
    for (rid, mis, terrace), q in refs.groupby(['rid', 'mis', 'terraza_id']):
        if len(q) < 2:
            continue
        cr = cruces_cota(perfiles[int(rid)], q.H_m.iloc[0])
        # No atribuir a B un salto de rama P_T: punto común aguas arriba de todos los B.
        cr = cr[(cr.s_PT_m > q.s_B_m.max()+1e-6) & cr.direccion.eq('ascendente')]
        if cr.empty:
            continue
        fixed = cr.iloc[0]
        for _, b in q.iterrows():
            rows.append({'rid': rid, 'mis': mis, 'terraza_id': terrace, 'caso_B': b.caso_B,
                'chi_B_m': b.chi_B_m, 's_B_m': b.s_B_m, 's_PT_comun_m': fixed.s_PT_m,
                'chi_PT_comun_m': fixed.chi_PT_m,
                'K_B_aislado_m02_anio': (fixed.chi_PT_m-b.chi_B_m)/(A0**THETA*b.T_ka*1000),
                'criterio': 'mismo H,T,P_T; primer ascendente común provisional'})
    return pd.DataFrame(rows)


def sensibilidad_altura(refs, perfiles):
    rows = []
    for _, b in refs.iterrows():
        for sign in [-1, 0, 1]:
            height = b.H_m+sign*b.error_H_publicado_m
            cr = cruces_cota(perfiles[int(b.rid)], height)
            cr = cr[(cr.s_PT_m > b.s_B_m+1e-6) & cr.direccion.eq('ascendente')]
            if cr.empty:
                rows.append({'caso_B': b.caso_B, 'cambio_H': sign, 'H_ensayo_m': height,
                             'estado': 'sin cruce ascendente'})
            for _, pt in cr.iterrows():
                rows.append({'caso_B': b.caso_B, 'cambio_H': sign, 'H_ensayo_m': height,
                    's_PT_m': pt.s_PT_m, 'chi_PT_m': pt.chi_PT_m,
                    'K_m02_anio': (pt.chi_PT_m-b.chi_B_m)/(A0**THETA*b.T_ka*1000),
                    'estado': 'escenario H ± magnitud publicada; no probabilístico'})
    return pd.DataFrame(rows)


def muestrear_swath(dem, linea, ancho_m=120., paso_m=12., paso_lateral_m=12.):
    if linea.geom_type != 'LineString' or linea.length <= 0 or ancho_m <= 0 or paso_m <= 0:
        raise ValueError('Transecto, ancho o paso inválido.')
    s = np.linspace(0, linea.length, max(2, int(np.ceil(linea.length/paso_m))+1))
    pts = np.array([linea.interpolate(float(a)).coords[0] for a in s])
    dxy = np.gradient(pts, axis=0)
    norm = np.linalg.norm(dxy, axis=1)
    if (norm == 0).any():
        raise ValueError('Transecto con inversión o nodos coincidentes; trazar una línea simple.')
    normal = np.column_stack([-dxy[:, 1], dxy[:, 0]])/norm[:, None]
    across = np.linspace(-ancho_m/2, ancho_m/2, max(3, int(np.ceil(ancho_m/paso_lateral_m))+1))
    xy = pts[:, None, :]+across[None, :, None]*normal[:, None, :]
    with rasterio.open(dem) as src:
        if src.crs != rasterio.crs.CRS.from_string(CRS):
            raise ValueError('DEM y capas deben estar en UTM18S; no se cambia la referencia vertical.')
        margin = max(src.res)*3
        xmin, ymin = xy.reshape(-1, 2).min(axis=0)-margin
        xmax, ymax = xy.reshape(-1, 2).max(axis=0)+margin
        win = from_bounds(xmin, ymin, xmax, ymax, src.transform).round_offsets().round_lengths()
        z = src.read(1, window=win, masked=True, boundless=True).astype(float).filled(np.nan)
        transform = src.window_transform(win)
        inv = ~transform
        cols = inv.a*xy[:, :, 0]+inv.b*xy[:, :, 1]+inv.c-.5
        rows = inv.d*xy[:, :, 0]+inv.e*xy[:, :, 1]+inv.f-.5
        values = map_coordinates(z, [rows, cols], order=1, mode='constant', cval=np.nan, prefilter=False)
    valid = np.isfinite(values).sum(axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        stats = np.nanpercentile(values, [0, 50, 90, 100], axis=1)
    # P90 es una envolvente exploratoria; no se denomina superficie original ni TerraceM.
    return pd.DataFrame({'distancia_m': s, 'x': pts[:, 0], 'y': pts[:, 1],
        'z_min': stats[0], 'z_mediana': stats[1], 'z_p90': stats[2], 'z_max': stats[3],
        'fraccion_cobertura': valid/values.shape[1], 'n_lateral': values.shape[1]})


def ajustar_quiebre(sw, transecto, curva='z_p90', cobertura_min=.9):
    intervals = ['plat_ini', 'plat_fin', 'acan_ini', 'acan_fin']
    vals = pd.to_numeric(pd.Series({k: transecto.get(k, np.nan) for k in intervals}), errors='coerce')
    if not np.isfinite(vals).all():
        return {'estado': 'pendiente: seleccionar intervalos de plataforma y paleoacantilado'}, None, None
    a, b, c, d = vals.to_numpy(float)
    if not (0 <= a < b < c < d <= sw.distancia_m.max()+1e-6):
        return {'estado': 'intervalos inválidos: orientar de mar a continente y separar plataforma/acantilado'}, None, None
    q = sw[np.isfinite(sw[curva]) & sw.fraccion_cobertura.ge(cobertura_min)]
    plat = q[q.distancia_m.between(a, b)]
    cliff = q[q.distancia_m.between(c, d)]
    if min(len(plat), len(cliff)) < 5:
        return {'estado': 'menos de cinco estaciones con cobertura suficiente en algún ajuste'}, None, None
    f1 = linregress(plat.distancia_m, plat[curva])
    f2 = linregress(cliff.distancia_m, cliff[curva])
    gap = f2.slope-f1.slope
    if gap <= 1e-8 or f2.slope <= 0:
        return {'estado': 'sin contraste plataforma-acantilado positivo; revisar orientación/selección'}, f1, f2
    x = (f1.intercept-f2.intercept)/gap
    H = f1.intercept+f1.slope*x
    # El quiebre debe quedar dentro de la extensión muestreada de ambos elementos.
    valid = a <= x <= d
    return {'estado': 'quiebre candidato, requiere revisión' if valid else 'intersección fuera de los elementos muestreados',
        's_quiebre_m': x, 'H_quiebre_DEM_m': H, 'pendiente_plataforma': f1.slope,
        'pendiente_acantilado': f2.slope, 'R2_plataforma': f1.rvalue**2,
        'R2_acantilado': f2.rvalue**2, 'curva': curva, 'aceptable_geometricamente': valid}, f1, f2


def revisar_transectos(capa, dem, output, curva='z_p90'):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    records = []
    if capa.empty:
        print('Aún no hay transectos: dibujarlos en QGIS, guardar ediciones y volver a esta celda.')
        return pd.DataFrame()
    if capa.id_perfil.duplicated().any():
        raise ValueError('id_perfil debe ser único en los transectos.')
    for _, tr in capa.iterrows():
        width = float(tr.ancho_m) if pd.notna(tr.ancho_m) else 120.
        sw = muestrear_swath(dem, tr.geometry, ancho_m=width)
        result, f1, f2 = ajustar_quiebre(sw, tr, curva=curva)
        fig, ax = plt.subplots(figsize=(10, 4.8), layout='constrained')
        ax.fill_between(sw.distancia_m, sw.z_min, sw.z_max, color='#dde3e7', label='Mínimo–máximo de banda')
        ax.plot(sw.distancia_m, sw.z_mediana, color='#64717b', label='Mediana')
        ax.plot(sw.distancia_m, sw.z_p90, color='#156c9c', label='Percentil 90 exploratorio')
        if f1 is not None:
            xjoint = result.get('s_quiebre_m', (tr.plat_fin+tr.acan_ini)/2)
            xx1 = np.array([tr.plat_ini, min(tr.acan_fin, max(tr.plat_fin, xjoint+25))])
            xx2 = np.array([max(tr.plat_ini, min(tr.acan_ini, xjoint-25)), tr.acan_fin])
            for fit, xx, col, lab in [(f1, xx1, '#19846e', 'Ajuste plataforma'), (f2, xx2, '#b05241', 'Ajuste paleoacantilado')]:
                ax.plot(xx, fit.intercept+fit.slope*xx, '--', color=col, label=lab)
            ax.axvspan(tr.plat_ini, tr.plat_fin, color='#19846e', alpha=.09)
            ax.axvspan(tr.acan_ini, tr.acan_fin, color='#b05241', alpha=.09)
            if 's_quiebre_m' in result and result.get('aceptable_geometricamente'):
                point = tr.geometry.interpolate(result['s_quiebre_m'])
                result.update(x_quiebre=point.x, y_quiebre=point.y)
                ax.scatter(result['s_quiebre_m'], result['H_quiebre_DEM_m'], color='black', s=35, zorder=5,
                           label='Intersección candidata')
        ax.set(xlabel='Distancia desde el inicio del transecto (m): mar → continente',
               ylabel='Elevación DEM (m)\nReferencia vertical por comprobar',
               title=f"{tr.id_perfil} · R{int(tr.rid):02d} · MIS {tr.mis} · margen {tr.margen}")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
        safe = re.sub(r'[^\w.-]+', '_', str(tr.id_perfil))
        fig.savefig(output/f'{safe}.png', dpi=180, facecolor='white')
        plt.show()
        plt.close(fig)
        sw.to_csv(output/f'{safe}_banda.csv', index=False)
        records.append({'id_perfil': tr.id_perfil, 'rid': tr.rid, 'mis': tr.mis, 'grupo': tr.grupo,
                        'margen': tr.margen, 'ancho_m': width} | result)
    return pd.DataFrame(records)


def fondo_relieve(ax, dem, bounds, max_lado=1200):
    xmin, ymin, xmax, ymax = bounds
    with rasterio.open(dem) as src:
        win = from_bounds(xmin, ymin, xmax, ymax, src.transform).round_offsets().round_lengths()
        h, w = int(win.height), int(win.width)
        factor = max(1., max(h, w)/max_lado)
        shape = (max(2, int(h/factor)), max(2, int(w/factor)))
        z = src.read(1, window=win, out_shape=shape, masked=True, boundless=True,
                     resampling=rasterio.enums.Resampling.bilinear).astype(float).filled(np.nan)
        if not np.isfinite(z).any():
            ax.set_facecolor('#e6edf2')
            return
        fill = np.where(np.isfinite(z), z, np.nanmedian(z))
        shade = LightSource(azdeg=315, altdeg=45).hillshade(fill, dx=(xmax-xmin)/shape[1], dy=(ymax-ymin)/shape[0])
        shade = np.ma.array(shade, mask=~np.isfinite(z))
        ax.imshow(shade, extent=[xmin, xmax, ymin, ymax], origin='upper', cmap='gray', vmin=0, vmax=1, alpha=.75)


def figuras_por_rio(refs, res, perfiles, nombres, capas, dem, output, margen_m=1200.):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for rid, refs_rio in refs.groupby('rid'):
        rid = int(rid)
        p = perfiles[rid]
        fig, axes = plt.subplots(2, 2, figsize=(16, 11), layout='constrained')
        ax = axes[0, 0]
        coords = [(float(v.x_B), float(v.y_B)) for _, v in refs_rio.iterrows()]
        for key in ['bordes', 'prolongaciones', 'transectos']:
            g = capas[key]
            if not g.empty:
                for line in g[g.rid.eq(rid)].geometry:
                    coords.extend(line.coords)
        xy = np.array(coords)
        lo, hi = xy.min(axis=0)-margen_m, xy.max(axis=0)+margen_m
        bounds = (lo[0], lo[1], hi[0], hi[1])
        fondo_relieve(ax, dem, bounds)
        ax.plot(p.x_utm, p.y_utm, color='#137aa2', lw=1.4, label='Cauce actual')
        if 'iniciales' in capas:
            g = capas['iniciales'].to_crs(CRS)
            g = g[pd.to_numeric(g['id']).eq(rid)]
            for j, (_, v) in enumerate(g.iterrows()):
                ax.plot(*v.geometry.xy, color='#94623e', ls=':', lw=1.1,
                        label='Línea inicial por cota' if j == 0 else None)
        if 'jm' in capas:
            jm = capas['jm']
            jm = jm[jm.x.between(lo[0], hi[0]) & jm.y.between(lo[1], hi[1])]
            for mis, marker in [('MIS 5e', '^'), ('MIS 5c', 's')]:
                pts = jm[jm.mis_asignado.eq(mis)]
                ax.scatter(pts.x, pts.y, marker=marker, facecolors='white', edgecolors='#44215f', s=30,
                           label=f'J-M {mis}')
        for key, col, ls, title in [('bordes', '#0e804e', '-', 'Borde observado (trazado propio)'),
                                     ('prolongaciones', '#bd2679', '--', 'Prolongación propia por el valle'),
                                     ('transectos', '#d18411', ':', 'Transecto de comprobación')]:
            g = capas[key]
            if not g.empty:
                g = g[g.rid.eq(rid)]
                for j, (_, v) in enumerate(g.iterrows()):
                    ax.plot(*v.geometry.xy, color=col, ls=ls, lw=1.7, label=title if j == 0 else None)
        colors = plt.get_cmap('tab10')
        subset = res[res.rid.eq(rid)] if not res.empty else pd.DataFrame()
        axchi, axdist, axk = axes[0, 1], axes[1, 0], axes[1, 1]
        axchi.plot(p.chi_m, p.elevacion_m, color='#697783', lw=1)
        axdist.plot(p.distancia_m/1000, p.elevacion_m, color='#697783', lw=1)
        ticklabels, yk = [], []
        for i, (_, b) in enumerate(refs_rio.iterrows()):
            col = colors(i % 10)
            lab = f"{b.mis}/{b.escenario}"
            ax.scatter(b.x_B, b.y_B, marker='x' if b.tipo_B == 'cota_inicial' else 'D', color=col, s=60,
                       label=f'B {lab}')
            ax.annotate(lab, (b.x_B, b.y_B), fontsize=7, xytext=(4, 4), textcoords='offset points')
            for a, xb in [(axchi, b.chi_B_m), (axdist, b.s_B_m/1000)]:
                a.scatter(xb, b.z_B_m, marker='X', color=col, s=40)
                a.axhline(b.H_m, color=col, alpha=.6, ls=':')
            q = subset[subset.caso_B.eq(b.caso_B)] if not subset.empty else pd.DataFrame()
            if not q.empty:
                for a, xcol, factor in [(axchi, 'chi_PT_m', 1), (axdist, 's_PT_m', 1000)]:
                    a.scatter(q[xcol]/factor, np.repeat(b.H_m, len(q)), facecolors='none', edgecolors=col, s=40)
                    for _, pt in q.iterrows():
                        a.annotate(str(int(pt.cruce_PT)), (pt[xcol]/factor, b.H_m), fontsize=7,
                                   xytext=(0, 5), textcoords='offset points')
                    primary = q[q.principal_operativo]
                    if not primary.empty:
                        u = primary.iloc[0]
                        a.scatter(u[xcol]/factor, b.H_m, color=col, s=25)
                        if a is axchi:
                            yarrow = b.H_m-max(5., abs(b.H_m-b.z_B_m)*.12)
                            a.annotate('', xy=(u.chi_PT_m, yarrow), xytext=(b.chi_B_m, yarrow),
                                       arrowprops={'arrowstyle':'<->','color':col,'lw':1})
                            a.text((u.chi_PT_m+b.chi_B_m)/2, yarrow,
                                   f'Δχ={u.delta_chi_m:.1f} m', color=col, fontsize=7, ha='center', va='bottom')
                valid = q[q.PT_elegible & q.K_George_m02_anio.gt(0)]
                if not valid.empty:
                    xx = np.full(len(valid), len(ticklabels), dtype=float)
                    xx += np.linspace(-.12, .12, len(valid))
                    axk.scatter(xx, valid.K_George_m02_anio, color=col, s=25)
                    axk.vlines(len(ticklabels), valid.K_George_m02_anio.min(), valid.K_George_m02_anio.max(), color=col)
                    ticklabels.append(lab)
                    yk.extend(valid.K_George_m02_anio)
        ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), xlabel='Este UTM18S (m)', ylabel='Norte UTM18S (m)',
               title='Reconstrucción local: puntos publicados y trazados propios')
        ax.set_aspect('equal')
        ax.ticklabel_format(useOffset=False, style='plain')
        ax.annotate('N', xy=(.95, .94), xytext=(.95, .79), xycoords='axes fraction', ha='center',
                    arrowprops={'arrowstyle': '-|>', 'color': 'black'})
        ax.legend(fontsize=6.7, loc='best', framealpha=.9)
        for a in [axchi, axdist]:
            a.grid(alpha=.2)
            a.set_ylabel('Elevación DEM (m)')
        axchi.set(xlabel='χ (m)', title='Todos los cruces de cota; número = alternativa P_T')
        xc = refs_rio.chi_B_m.to_list()
        if not subset.empty:
            xc += subset.chi_PT_m.to_list()
        xlo, xhi = min(xc), max(xc)
        margin = max(25., (xhi-xlo)*.15)
        axchi.set_xlim(max(float(p.chi_m.iloc[0]), xlo-margin), min(float(p.chi_m.iloc[-1]), xhi+margin))
        axchi.set_ylim(min(float(refs_rio.z_B_m.min()), float(refs_rio.H_m.min()))-10,
                      float(refs_rio.H_m.max())+max(10, float(refs_rio.H_m.max()-refs_rio.z_B_m.min())*.15))
        axdist.set(xlabel='Distancia aguas arriba (km)', title='Perfil longitudinal completo')
        axk.set_xticks(range(len(ticklabels)), ticklabels, rotation=30, ha='right', fontsize=7)
        axk.set(ylabel='K condicional (m⁰·²/año)', title='Variación entre P_T con cada B fijo; línea = rango')
        axk.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
        axk.grid(alpha=.2, axis='y')
        if not yk:
            axk.text(.5, .5, 'Sin cruce ascendente utilizable', transform=axk.transAxes, ha='center')
        fig.suptitle(f'R{rid:02d} {nombres[rid]} · geometría B–PT · H DEM de la terraza · todos los escenarios', fontsize=13)
        fig.savefig(output/f'R{rid:02d}_B_PT_K.png', dpi=190, facecolor='white')
        plt.close(fig)


def exportar_corrida(folder, cfg, tablas, fuentes):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for name, table in tablas.items():
        table.to_csv(folder/f'{name}.csv', index=False, encoding='utf-8-sig')
    hashes = {}
    for path in fuentes:
        p = Path(path)
        if p.is_file():
            hashes[str(p.resolve())] = huella(p)
    # Los shapefiles requieren DBF/SHX/PRJ para conservar geometría y atributos.
    for path in list(fuentes):
        p = Path(path)
        if p.suffix.lower() == '.shp':
            for ext in ['.shx', '.dbf', '.prj', '.cpg']:
                side = p.with_suffix(ext)
                if side.is_file():
                    hashes[str(side.resolve())] = huella(side)
    (folder/'Trazabilidad.json').write_text(json.dumps({'fecha_utc': datetime.now(timezone.utc).isoformat(),
        'config': cfg, 'fuentes_sha256': hashes,
        'alcance': 'calibración condicional; sin validación automática de B, P_T ni cronología',
        'errores': 'magnitudes publicadas; aproximación local, sin B ni error DEM; no IC total'},
        ensure_ascii=False, indent=2, default=str), encoding='utf8')
    return folder


def copiar_gpkg(origen, destino):
    """Snapshot coherente que incluye transacciones confirmadas en WAL de SQLite."""
    import sqlite3
    with sqlite3.connect(Path(origen).resolve().as_uri()+'?mode=ro', uri=True) as source:
        with sqlite3.connect(str(destino)) as target:
            source.backup(target)
