"""Funciones compartidas. Entrega estática: no se ha ejecutado el análisis."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import unicodedata
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.stats import linregress

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
CFG_PATH = BASE / 'Configuracion.json'
CFG = json.loads(CFG_PATH.read_text(encoding='utf-8'))
RUN_ID = CFG['run_id']
if not re.fullmatch(r'v[1-9]\d*', RUN_ID):
    raise ValueError('run_id debe ser v1, v2, etc.')
EXEC_DATE = RUN_ID
for key in ['dem','kp_manual_csv','terrazas_excel','revision_muestras_excel',
            'seleccion_excel','geologia_integrada','fallas_integradas','chaf','maule',
            'outlets_csv','originales_18_csv','fuente_JM_csv','escenarios_30_csv',
            'grant_notebook','graben_shp','topoanalysis_codigo']:
    CFG[key] = str((BASE / CFG[key]).resolve())
CFG['base_dir'] = str(BASE)
ROOT = BASE / 'Productos' / 'Datos' / RUN_ID
_pathmap = BASE/'Archivos_Memoria/Rutas_fuentes.json'
PATHMAP = json.loads(_pathmap.read_text(encoding='utf-8')) if _pathmap.exists() else {}
def fuente(path):
    if Path(path).exists():
        return Path(path)
    key = str(path).replace('\\','/').lower()
    if key not in PATHMAP:
        raise FileNotFoundError('Dependencia no incluida: ' + str(path))
    return BASE / PATHMAP[key]

def version_path(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r'_v[1-9]\d*$', '', path.stem)
    v = 1
    while (path.parent / f'{stem}_v{v}{path.suffix}').exists(): v += 1
    return path.parent / f'{stem}_v{v}{path.suffix}'

def nueva_carpeta(folder):
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=True)
    v = 1
    while (folder / f'v{v}').exists(): v += 1
    result = folder / f'v{v}'; result.mkdir()
    return result

def guardar_figura(fig, tipo, rid=None, carpeta='02_Perfiles'):
    prefix = f'{int(rid):02d}_' if rid is not None else ''
    path = version_path(BASE/'Productos'/carpeta/(prefix+tipo+'.png'))
    while path.with_suffix('.pdf').exists():
        stem=path.stem.rsplit('_v',1)[0]
        v=int(path.stem.rsplit('_v',1)[1])+1
        path=path.with_name(f'{stem}_v{v}{path.suffix}')
    fig.savefig(path, dpi=300, bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'), bbox_inches='tight')
    return path

EXPECTED = [str(i) for i in range(1, 19)]
COLS = ['distancia_m', 'elevacion_m', 'area_m2', 'x_utm', 'y_utm']

def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def signature(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()

def input_hash(path):
    path = Path(path)
    if path.suffix.lower() == '.shp':
        return signature({p.name: file_hash(p) for p in sorted(path.parent.glob(path.stem + '.*'))})
    return file_hash(path)

def profile_hash(p):
    return hashlib.sha256(np.asarray(p[COLS].to_numpy(float), dtype='<f8', order='C').tobytes()).hexdigest()

def normalize_id(value):
    number = float(value)
    if not np.isfinite(number) or not number.is_integer() or str(int(number)) not in EXPECTED:
        raise ValueError(f'ID de río inválido: {value}')
    return str(int(number))

def text(value):
    return '' if pd.isna(value) else str(value).strip()

def normalizar_nombre(value):
    value = ''.join(c for c in unicodedata.normalize('NFKD', text(value)) if not unicodedata.combining(c))
    value = re.sub(r'^s/n\s*', '', value.lower().strip())
    value = re.sub(r'^(rio|estero)\s+', '', value)
    value = re.sub(r'\s+', ' ', value)
    return {'buchupero': 'buchupureo', 'cobequecura': 'cobquecura'}.get(value, value)

def output(stage, name, extension='csv'):
    folder = ROOT / stage
    folder.mkdir(parents=True, exist_ok=True)
    # El manifiesto es el índice explícito de la corrida. Los productos nuevos
    # llevan además la fecha real de esta sesión, sin buscar el "último" archivo.
    stamp = RUN_ID
    return folder / f'{name}_{stamp}.{extension}'

def save_table(stage, name, table):
    path = output(stage, name)
    table.to_csv(path, index=False, encoding='utf-8-sig')
    return path

def save_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(data, ensure_ascii=False, indent=2, default=_json), encoding='utf-8')
    temp.replace(path)

def _json(value):
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, Path): return str(value)
    raise TypeError(type(value).__name__)

def begin(stage):
    """Un intento interrumpido no debe dejar un manifiesto marcado completo."""
    save_json(output(stage, 'manifest', 'json'), {'run_id': RUN_ID, 'estado': 'en_proceso'})

def finish(stage, files, dependencies=None, extra=None):
    files = {k: {'path': str(Path(v).resolve()), 'sha256': file_hash(v)} for k, v in files.items()}
    notebook = globals().get('NOTEBOOK_PATH')
    document = {'run_id': RUN_ID, 'estado': 'completo', 'fecha_utc': datetime.now(timezone.utc).isoformat(),
                'files': files, 'dependencies': dependencies or {}, 'config': CFG,
                'codigo_sha256': file_hash(__file__), 'notebook': str(notebook) if notebook else None,
                'notebook_sha256': file_hash(notebook) if notebook and Path(notebook).is_file() else None,
                **(extra or {})}
    save_json(output(stage, 'manifest', 'json'), document)
    return document

def read_stage(stage):
    path = output(stage, 'manifest', 'json')
    if not path.is_file():
        raise FileNotFoundError(f'Ejecutar primero {stage}: {path}')
    m = json.loads(path.read_text(encoding='utf-8'))
    if m.get('estado') != 'completo' or m.get('run_id') != RUN_ID:
        raise ValueError(f'{stage}: corrida incompleta o incompatible.')
    # Cada dependencia se referencia por archivo y huella, nunca por "último".
    for filename, expected in m.get('dependencies', {}).items():
        if file_hash(filename) != expected:
            raise ValueError(f'{stage}: entrada modificada; recalcular: {filename}')
        dep = Path(filename)
        if dep.name == f'manifest_{RUN_ID}.json' and dep.parent.parent == ROOT:
            read_stage(dep.parent.name)
    for item in m['files'].values():
        if file_hash(item['path']) != item['sha256']:
            raise ValueError(f'{stage}: producto modificado: {item["path"]}')
    return m

def dependency(stage):
    read_stage(stage)
    path = output(stage, 'manifest', 'json')
    return {str(path): file_hash(path)}

def stage_table(stage, key):
    m = read_stage(stage)
    return pd.read_csv(m['files'][key]['path'], dtype={'ID': str, 'rio_id': str})

def load_profiles():
    m = read_stage('perfiles')
    for key in ['theta', 'A0_m2', 'area_minima_m2', 'ventana_m', 'min_puntos', 'metodo_chi']:
        if m['config'].get(key) != CFG[key]:
            raise ValueError(f'Configuración de perfiles modificada ({key}); recalcular o restaurar configuración.')
    if set(m['rios']) != set(EXPECTED):
        raise ValueError('La exportación de perfiles no contiene los 18 IDs.')
    profiles, names = {}, {}
    for rid, item in m['rios'].items():
        p = pd.read_csv(m['files'][f'rio_{rid}']['path'], float_precision='round_trip')
        if set(COLS) - set(p): raise ValueError(f'Río {rid}: faltan columnas.')
        s, a = p['distancia_m'].to_numpy(float), p['area_m2'].to_numpy(float)
        if len(p) < 2 or not np.isfinite(s).all() or not np.isclose(s[0], 0) or np.any(np.diff(s) <= 0):
            raise ValueError(f'Río {rid}: distancia inválida.')
        if not np.isfinite(a).all() or np.any(a <= 0): raise ValueError(f'Río {rid}: área inválida.')
        if not np.isfinite(p[['x_utm', 'y_utm']].to_numpy(float)).all(): raise ValueError(f'Río {rid}: coordenadas inválidas.')
        p.attrs.update(item['parametros'])
        if profile_hash(p) != item['perfil_hash']: raise ValueError(f'Río {rid}: huella incompatible.')
        profiles[rid], names[rid] = p, item['nombre']
    return profiles, names, m

def chi_1d(p, theta=None, A0=None):
    theta = CFG['theta'] if theta is None else float(theta)
    A0 = CFG['A0_m2'] if A0 is None else float(A0)
    s, a = p['distancia_m'].to_numpy(float), p['area_m2'].to_numpy(float)
    if not np.isfinite([theta, A0]).all() or A0 <= 0 or theta < 0:
        raise ValueError('Parámetros de chi inválidos.')
    if not np.isfinite(s).all() or not np.isfinite(a).all() or np.any(np.diff(s) <= 0) or np.any(a <= 0):
        raise ValueError('Chi requiere distancias crecientes y áreas positivas, sin huecos.')
    return cumulative_trapezoid((A0 / a)**theta, x=s, initial=0.0)

def slope_ksn(p, window, minimum):
    """No ajustar a través de huecos de elevación ni usar abs(S)."""
    if not np.isfinite(window) or window <= 0 or minimum < 3: raise ValueError('Ventana inválida.')
    result = p.copy()
    s, z, a = (p[c].to_numpy(float) for c in ['distancia_m', 'elevacion_m', 'area_m2'])
    slope, r2 = np.full(len(p), np.nan), np.full(len(p), np.nan)
    count = np.zeros(len(p), dtype=int)
    for j, center in enumerate(s):
        lo, hi = np.searchsorted(s, [center - window / 2, center + window / 2])
        hi = np.searchsorted(s, center + window / 2, side='right')
        valid = np.isfinite(z[lo:hi])
        count[j] = int(valid.sum())
        if valid.all() and count[j] >= minimum:
            fit = linregress(s[lo:hi] - center, z[lo:hi])
            slope[j], r2[j] = fit.slope, fit.rvalue**2
    good = np.isfinite(slope) & (slope > 0) & np.isfinite(a) & (a > 0)
    result['pendiente'], result['r2_pendiente'] = slope, r2
    result['ksn'] = np.where(good, slope * a**CFG['theta'], np.nan)
    result['ksn_valido'] = good
    result['n_puntos_ventana'] = count
    result['ventana_parcial'] = (s - window/2 < s[0]) | (s + window/2 > s[-1])
    result['pendiente_no_positiva'] = np.isfinite(slope) & (slope <= 0)
    result.attrs.update(p.attrs, ventana_m=window, min_puntos=minimum, politica_huecos='rechazar_ventana')
    return result

def load_manuals(profiles):
    raw = pd.read_csv(CFG['kp_manual_csv'], sep=CFG['kp_separador'], decimal=CFG['kp_decimal'], dtype={'ID':str}).dropna(how='all')
    raw.columns = raw.columns.str.strip()
    good, bad = [], []
    seen = set()
    for r in raw.to_dict('records'):
        try:
            rid=normalize_id(r['ID'] if pd.notna(r['ID']) else r.get('id')); r['ID']=rid
            if not text(r['id_kp']):
                if any(pd.notna(r.get(k)) for k in ['x_utm','y_utm','elevacion_m']):
                    raise ValueError('Coordenadas sin id_kp')
                r.update(estado='sin_seleccion',detalle='sin KP manual'); bad.append(r); continue
            if r['id_kp'] in seen: raise ValueError('KP duplicado')
            seen.add(r['id_kp']); p=profiles[rid]
            xy=np.array([r['x_utm'],r['y_utm']],float)
            dist=np.linalg.norm(p[['x_utm','y_utm']].to_numpy(float)-xy,axis=1)
            ids=np.flatnonzero(dist<=0.10)
            if len(ids)!=1: raise ValueError('No hay un nodo XY único dentro de 0.10 m')
            j=int(ids[0]); dz=abs(float(p.elevacion_m.iloc[j])-float(r['elevacion_m']))
            if not np.isfinite(dz) or dz>0.10: raise ValueError('Elevación incompatible')
            original=r['perfil_hash']; h=profile_hash(p)
            r.update(indice_original=r['indice'],distancia_original_m=r['distancia_m'],
                hash_original=original,perfil_hash=h,indice=j,distancia_m=float(p.distancia_m.iloc[j]),
                elevacion_m=float(p.elevacion_m.iloc[j]),desajuste_xy_m=float(dist[j]),desajuste_z_m=dz,
                metodo='Manual',estado='ok' if original==h else 'manual_reubicado',detalle='')
            good.append(r)
        except (ValueError,KeyError,TypeError) as exc:
            r.update(estado='incompatible',detalle=str(exc));bad.append(r)
    cols=list(dict.fromkeys([*raw.columns,'estado','detalle','metodo','hash_original','indice_original',
        'distancia_original_m','desajuste_xy_m','desajuste_z_m']))
    accepted=pd.DataFrame(good,columns=cols); rejected=pd.DataFrame(bad,columns=cols)
    if not rejected.empty and rejected.estado.eq('incompatible').any():
        raise ValueError('KP incompatibles; revisar CSV y geometría: '+rejected.loc[rejected.estado.eq('incompatible'),['ID','detalle']].to_string(index=False))
    return accepted,rejected


def crossings(p, level):
    s, z = p['distancia_m'].to_numpy(float), p['elevacion_m'].to_numpy(float)
    points, flats = [], []
    for j in range(len(s)-1):
        if not np.isfinite(z[j:j+2]).all(): continue
        if z[j] == z[j+1]:
            if z[j] == level: flats.append((s[j], s[j+1]))
        elif min(z[j], z[j+1]) <= level <= max(z[j], z[j+1]):
            points.append(s[j] + (level-z[j])/(z[j+1]-z[j])*(s[j+1]-s[j]))
    return np.unique(np.round(points, 6)), flats

def interval_fit(p, start, stop, theta, A0):
    s, z = p['distancia_m'].to_numpy(float), p['elevacion_m'].to_numpy(float)
    if not np.isfinite([start, stop]).all() or not s[0] <= start < stop <= s[-1]: raise ValueError('Límites fuera del perfil')
    # Incluir los nodos que soportan los extremos para no interpolar sobre huecos.
    lo = max(0, np.searchsorted(s, start, side='right')-1)
    hi = min(len(s)-1, np.searchsorted(s, stop))
    if not np.isfinite(z[lo:hi+1]).all(): raise ValueError('Intervalo con huecos de elevación')
    idx = (s >= start) & (s <= stop)
    if idx.sum() < 10 or stop-start < 500: raise ValueError('Calibración requiere 10 nodos y 500 m como mínimo')
    chi = chi_1d(p, theta, A0)
    fit = linregress(chi[idx], z[idx])
    if not np.isfinite(fit.slope) or fit.slope <= 0: raise ValueError('Pendiente de calibración no positiva')
    return {'ksn': float(fit.slope*A0**theta), 'R2': float(fit.rvalue**2), 'n_puntos': int(idx.sum()), 'chi': chi}
