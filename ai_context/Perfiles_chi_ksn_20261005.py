# Code-only reading snapshot; not an independently executable script.
# Notebook is authoritative. Cell magics, external data and cell order still apply.

# Notebook cell 2
from pathlib import Path
import sys
import re
CARPETA_NOTEBOOKS = Path.cwd()
if not (CARPETA_NOTEBOOKS / 'post_auditoria_20261005.py').is_file():
    CARPETA_NOTEBOOKS = Path(r'C:/datos/Memoria/Codigos Stanford/Notebook_post_auditoria')
sys.path.insert(0, str(CARPETA_NOTEBOOKS))
import post_auditoria_20261005 as pa
import importlib
pa = importlib.reload(pa)
from post_auditoria_20261005 import CFG, RUN_ID
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
print('Corrida compartida:', RUN_ID)

pa.NOTEBOOK_PATH = CARPETA_NOTEBOOKS / 'Perfiles_chi_ksn_20261005.ipynb'


# Notebook cell 4
from TopoAnalysis import dem as d
import rasterio
from rasterio.enums import Resampling
from scipy.stats import linregress
topo = pa.read_stage('topografia')
deps = pa.dependency('topografia')
DEM_PATH = Path(CFG['dem'])
if pa.file_hash(DEM_PATH) != topo['dem_sha256']: raise ValueError('El DEM cambió desde topografía.')
deps[str(DEM_PATH)] = topo['dem_sha256']
AREA_PATH = Path(topo['files']['area']['path'])
FD_PATH = Path(topo['files']['flow_direction']['path'])
catalogo = pd.read_csv(topo['files']['outlets']['path'], dtype={'rio_id': str})
river_ids = dict(zip(catalogo.rio_id, catalogo.nombre))
sample_outlets = {r.nombre: (r.x, r.y) for r in catalogo.itertuples(index=False)}
A0, AREA_THRESHOLD, THETA_REF = CFG['A0_m2'], CFG['area_minima_m2'], CFG['theta']
colors = {name: plt.get_cmap('tab20')(i) for i, name in enumerate(sample_outlets)}
pa.begin('perfiles')


# Notebook cell 5
# Comprobar que existan las entradas
for path in (AREA_PATH, FD_PATH, DEM_PATH):
    if not path.exists():
        raise FileNotFoundError(f"No se encontró: {path}")


# ── Cargar las grillas con TopoAnalysis ───────────────────────────────────
area = d.Area.load(str(AREA_PATH))
fd = d.FlowDirectionD8.load(str(FD_PATH))
elevation = d.Elevation.load(str(DEM_PATH))


# ── Acceder a las matrices ───────────────────────────────────────────────
area_grid = area._griddata
fd_grid = fd._griddata
elev_grid = elevation._griddata


# ── Comprobar dimensiones ────────────────────────────────────────────────
if not (area_grid.shape == fd_grid.shape == elev_grid.shape):
    raise ValueError(
        "Las grillas tienen dimensiones diferentes:\n"
        f"Área: {area_grid.shape}\n"
        f"Dirección de flujo: {fd_grid.shape}\n"
        f"Elevación: {elev_grid.shape}"
    )

print("Grillas cargadas correctamente.")
print(f"Dimensiones: {elev_grid.shape[0]} filas × "
      f"{elev_grid.shape[1]} columnas")
print(f"DEM de elevación: {DEM_PATH.name}")
print("Siguiente paso: comprobar la alineación espacial y los outlets.")


# Notebook cell 6
import rasterio

# ── Leer la información espacial ─────────────────────────────────────────
grid_info = {}

for name, path in {
    "Elevación": DEM_PATH,
    "Área": AREA_PATH,
    "Dirección": FD_PATH,
}.items():
    with rasterio.open(path) as src:
        grid_info[name] = {
            "crs": src.crs,
            "transform": src.transform,
            "shape": src.shape,
            "res": src.res,
        }

reference = grid_info["Elevación"]

if reference["crs"] is None:
    raise ValueError("El DEM no tiene un sistema de coordenadas definido.")

if not reference["crs"].is_projected:
    raise ValueError("El DEM debe estar en coordenadas proyectadas.")

for name, info in grid_info.items():
    if info["crs"] != reference["crs"]:
        raise ValueError(f"{name}: sistema de coordenadas diferente.")

    if info["shape"] != reference["shape"]:
        raise ValueError(f"{name}: dimensiones diferentes.")

    if not info["transform"].almost_equals(reference["transform"]):
        raise ValueError(f"{name}: alineación o resolución diferente.")

print("Las tres grillas coinciden espacialmente.")
print(f"Sistema de coordenadas: {reference['crs']}")
print(f"Resolución: {reference['res']}")
print(f"Unidad: {reference['crs'].linear_units}")


# Notebook cell 8
# ── Función para extraer un perfil longitudinal ──────────────────────────
def extraer_perfil(outlet, area_minima=AREA_THRESHOLD):
    row, col = fd._xy_to_rowscols((outlet,))[0]

    if row is None or col is None:
        raise ValueError("El outlet está fuera de la grilla.")

    nrows, ncols = area_grid.shape
    recorrido = []
    visitadas = set()

    while True:
        if (row, col) in visitadas:
            raise ValueError("Se detectó un ciclo en la dirección de flujo.")

        visitadas.add((row, col))
        drainage_area = float(area_grid[row, col])

        if not np.isfinite(drainage_area) or drainage_area < area_minima:
            if not recorrido:
                raise ValueError(
                    "El outlet no tiene área válida o está bajo el umbral."
                )
            break

        recorrido.append((row, col))

        # Vecinos que realmente drenan hacia la celda actual
        candidatos = []

        for r, c in fd.get_upstream_cell_indexes(row, col):
            r, c = int(r), int(c)

            if 0 <= r < nrows and 0 <= c < ncols:
                a = float(area_grid[r, c])

                if np.isfinite(a) and a >= area_minima:
                    candidatos.append((r, c))

        if not candidatos:
            break

        # Criterio de cauce principal: mayor área aportante
        row, col = max(
            candidatos,
            key=lambda rc: float(area_grid[rc[0], rc[1]])
        )

    if len(recorrido) < 2:
        raise ValueError(
            "El perfil tiene menos de dos puntos. Revisa el outlet y el umbral."
        )

    rows, cols = np.array(recorrido, dtype=int).T
    xy = np.asarray(fd._rowscols_to_xy(recorrido), dtype=float)
    x, y = xy[:, 0], xy[:, 1]

    # Distancia real entre centros de celdas, incluidas las diagonales
    pasos = np.hypot(np.diff(x), np.diff(y))
    distancia = np.concatenate(([0.0], np.cumsum(pasos)))

    perfil = pd.DataFrame({
        "fila": rows,
        "columna": cols,
        "x_utm": x,
        "y_utm": y,
        "distancia_m": distancia,
        "elevacion_m": elev_grid[rows, cols],
        "area_m2": area_grid[rows, cols],
    })

    # Conservar las posiciones sin elevación para no unir huecos en el gráfico
    perfil["elevacion_valida"] = np.isfinite(perfil["elevacion_m"])
    perfil.loc[~perfil["elevacion_valida"], "elevacion_m"] = np.nan

    # Diferencia entre la coordenada indicada y el centro de su celda
    perfil.attrs["ajuste_a_celda_m"] = float(
        np.hypot(x[0] - outlet[0], y[0] - outlet[1])
    )

    return perfil


# Notebook cell 9
perfiles_ksn, sensibilidad_ksn, estados = {}, {}, []
archivos, rios = {}, {}
for rid in pa.EXPECTED:
    registro = {'ID': rid, 'rio': river_ids[rid], 'estado': 'error', 'detalle': ''}
    try:
        p = extraer_perfil(sample_outlets[river_ids[rid]], area_minima=AREA_THRESHOLD)
        p['chi_m'] = pa.chi_1d(p)
        # Diagnóstico de discretización, sin afirmar que reproduce d.Chi.
        s, a = p['distancia_m'].to_numpy(float), p['area_m2'].to_numpy(float)
        rect = np.r_[0., np.cumsum(np.diff(s)*(A0/a[:-1])**THETA_REF)]
        p.attrs.update(theta_ref=THETA_REF, A0_m2=A0, area_minima_m2=AREA_THRESHOLD,
                       metodo_chi=CFG['metodo_chi'], run_id=RUN_ID)
        base = pa.slope_ksn(p, CFG['ventana_m'], CFG['min_puntos'])
        perfiles_ksn[rid] = base
        for window in sorted(set([CFG['ventana_m'], *CFG['ventanas_sensibilidad_m']])):
            result = base if window == CFG['ventana_m'] else pa.slope_ksn(p, window, CFG['min_puntos'])
            sensibilidad_ksn[(rid, window)] = result
            archivos[f'sens_{rid}_{window:g}'] = pa.save_table('perfiles', f'sens_rio_{int(rid):02d}_{window:g}m', result)
        archivos[f'rio_{rid}'] = pa.save_table('perfiles', f'rio_{int(rid):02d}', base)
        rios[rid] = {'nombre': river_ids[rid], 'archivo': archivos[f'rio_{rid}'].name,
                     'parametros': dict(base.attrs), 'perfil_hash': pa.profile_hash(base)}
        registro.update(estado='ok', puntos=len(p), sin_elevacion=int(p.elevacion_m.isna().sum()),
            ksn_validos=int(base.ksn_valido.sum()), pendientes_no_positivas=int(base.pendiente_no_positiva.sum()),
            ventanas_parciales=int(base.ventana_parcial.sum()),
            diferencia_trapecio_rectangulo_aval_m=float(np.max(np.abs(p.chi_m.to_numpy()-rect))))
    except (ValueError, IndexError, RuntimeError) as exc:
        registro['detalle'] = str(exc)
    estados.append(registro)
    archivos['qc'] = pa.save_table('perfiles', 'qc_perfiles', pd.DataFrame(estados))
resumen_qc = pd.DataFrame(estados)
display(resumen_qc)
if set(perfiles_ksn) != set(pa.EXPECTED):
    raise RuntimeError('Faltan perfiles: revisar qc_perfiles antes de continuar.')


# Notebook cell 11
import geopandas as gpd
from shapely.geometry import LineString
cauces = gpd.GeoDataFrame([{'ID': rid, 'rio': river_ids[rid], 'perfil_hash': pa.profile_hash(p),
    'geometry': LineString(p[['x_utm','y_utm']].to_numpy(float))} for rid,p in perfiles_ksn.items()], crs=CFG['crs'])
ruta_cauces = pa.output('perfiles', 'cauces_principales', 'gpkg')
cauces.to_file(ruta_cauces, layer='cauces', driver='GPKG', mode='w')
archivos['cauces_principales'] = ruta_cauces
pa.finish('perfiles', archivos, deps, extra={'rios': rios, 'unidad_distancia': 'm', 'unidad_area': 'm2',
    'metodo_chi': CFG['metodo_chi'], 'criterio_cauce': 'Mayor área aportante en cada confluencia'})


# Notebook cell 13
RIOS_VER = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]
GUARDAR_FIGURAS = False
for rid in map(str, RIOS_VER):
    p = perfiles_ksn[rid]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].plot(p.distancia_m/1000, p.elevacion_m)
    ax[1].plot(p.chi_m, p.elevacion_m)
    for (i,w), result in sensibilidad_ksn.items():
        if i == rid: ax[2].plot(result.distancia_m/1000, result.ksn, label=f'{w:g} m')
    ax[0].set(xlabel='Distancia (km)', ylabel='Elevación (m)')
    ax[1].set(xlabel='χ (m)', ylabel='Elevación (m)')
    ax[2].set(xlabel='Distancia (km)', ylabel=f'ksn (m^{2*THETA_REF:g})')
    ax[2].legend(); fig.suptitle(river_ids[rid]); fig.tight_layout()
    if GUARDAR_FIGURAS: fig.savefig(pa.output('perfiles', f'fig_rio_{int(rid):02d}', 'png'), dpi=300)
    plt.show(); plt.close(fig)


# Notebook cell 14

