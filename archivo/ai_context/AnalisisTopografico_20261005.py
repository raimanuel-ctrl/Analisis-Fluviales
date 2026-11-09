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

pa.NOTEBOOK_PATH = CARPETA_NOTEBOOKS / 'AnalisisTopografico_20261005.ipynb'


# Notebook cell 4
from TopoAnalysis import dem as d
from pathlib import Path
from datetime import datetime
import json
import sys
import importlib.metadata as metadata

import numpy as np
import pandas as pd
import rasterio
from rasterio.features import shapes
from rasterio.windows import Window
import geopandas as gpd
from shapely.geometry import shape

BASE_DIR = Path(CFG['base_dir'])
DEM_PATH = Path(CFG['dem'])
FECHA = RUN_ID + '_' + pa.EXEC_DATE
DATA_DIR = pa.ROOT / 'topografia'
MASK_DIR = DATA_DIR / 'masks'
POLY_DIR = DATA_DIR / 'basin_polygons'
CRS_OUTLETS = rasterio.crs.CRS.from_string(CFG['crs'])
GUARDAR_LOGAREA = False
REUTILIZAR_HIDROLOGIA = True


# Notebook cell 5
catalogo = pd.read_csv(CARPETA_NOTEBOOKS / 'outlets_20261005.csv')
catalogo['rio_id'] = catalogo['rio_id'].map(pa.normalize_id)
if len(catalogo) != 18 or set(catalogo['rio_id']) != set(pa.EXPECTED):
    raise ValueError('Se requieren los 18 IDs, sin duplicados.')
if catalogo['nombre'].duplicated().any() or catalogo[['x', 'y']].duplicated().any():
    raise ValueError('Catálogo de outlets duplicado.')
if not np.isfinite(catalogo[['x', 'y']].to_numpy(float)).all():
    raise ValueError('Coordenadas inválidas.')
if not catalogo['crs'].eq(CFG['crs']).all(): raise ValueError('CRS del catálogo incompatible.')
samples = catalogo[['nombre', 'x', 'y']].values.tolist()
display(catalogo)


# Notebook cell 7
if not DEM_PATH.is_file():
    raise FileNotFoundError(f'No existe el DEM: {DEM_PATH}')

controles = []
with rasterio.open(DEM_PATH) as src:
    if src.count != 1:
        raise ValueError('Se esperaba un DEM de una sola banda.')
    if src.crs is None or not src.crs.is_projected:
        raise ValueError('El DEM debe tener un CRS proyectado definido.')
    if not np.isclose(src.crs.linear_units_factor[1], 1.0):
        raise ValueError('Las unidades horizontales deben ser metros.')
    if src.crs != CRS_OUTLETS:
        raise ValueError(f'CRS incompatibles: DEM={src.crs}; outlets={CRS_OUTLETS}.')

    transform_dem = src.transform
    if (transform_dem.b != 0 or transform_dem.d != 0
            or transform_dem.a <= 0 or transform_dem.e >= 0
            or not np.isclose(transform_dem.a, -transform_dem.e, rtol=1e-9, atol=1e-9)):
        raise ValueError('Se requiere una grilla al norte, sin rotación y con píxeles cuadrados.')

    crs_dem = src.crs
    shape_dem = src.shape
    resolucion_m = float(transform_dem.a)
    area_celda_m2 = abs(transform_dem.a * transform_dem.e - transform_dem.b * transform_dem.d)
    n_total = src.width * src.height
    n_validas = 0
    for _, ventana in src.block_windows(1):
        bloque = src.read(1, window=ventana, masked=True)
        validas = ~np.ma.getmaskarray(bloque) & np.isfinite(bloque.data)
        n_validas += int(np.count_nonzero(validas))
    del bloque, validas
    if n_validas == 0:
        raise ValueError('El DEM no contiene elevaciones válidas.')

    xllcenter = transform_dem.c + resolucion_m / 2
    yllcenter = transform_dem.f - resolucion_m * (src.height - 0.5)
    for rio in catalogo.itertuples(index=False):
        columna = int(round((rio.x - xllcenter) / resolucion_m))
        fila = int((src.height - 1) - round((rio.y - yllcenter) / resolucion_m))
        detalle = ''
        if not (0 <= fila < src.height and 0 <= columna < src.width):
            detalle = 'Outlet fuera de la grilla'
        else:
            valor = src.read(1, window=Window(columna, fila, 1, 1), masked=True)
            if valor.count() != 1 or not np.isfinite(valor.compressed()).all():
                detalle = 'Outlet sin elevación válida'
        controles.append({'rio_id': rio.rio_id, 'nombre': rio.nombre,
                          'fila': fila, 'columna': columna,
                          'estado': 'error' if detalle else 'ok', 'detalle': detalle})

for carpeta in (DATA_DIR, MASK_DIR, POLY_DIR):
    carpeta.mkdir(parents=True, exist_ok=True)
catalogo.to_csv(DATA_DIR / f'catalogo_outlets_{FECHA}.csv', index=False, encoding='utf-8-sig')
qc_outlets = pd.DataFrame(controles)
qc_outlets.to_csv(DATA_DIR / f'qc_outlets_{FECHA}.csv', index=False, encoding='utf-8-sig')
if not qc_outlets['estado'].eq('ok').all():
    display(qc_outlets)
    raise ValueError('Hay outlets inválidos. Consultar qc_outlets antes de continuar.')

print('CRS:', crs_dem)
print(f'Resolución: {resolucion_m:g} m')
print(f'Celdas válidas: {n_validas:,} / {n_total:,}')
print(f'Sin datos: {100 * (1 - n_validas / n_total):.2f}%')
print('Los 18 outlets tienen cobertura válida.')


# Notebook cell 9
from importlib import metadata
codigo_topo = Path(d.__file__).parent
implementacion = {str(p.relative_to(codigo_topo)): pa.file_hash(p)
                  for p in sorted(codigo_topo.rglob('*'))
                  if p.is_file() and p.suffix.lower() in {'.py', '.pyd', '.so'}}
versiones_cache = {}
for paquete in ['numpy', 'scipy', 'GDAL', 'TopoAnalysis']:
    try: versiones_cache[paquete] = metadata.version(paquete)
    except metadata.PackageNotFoundError: versiones_cache[paquete] = 'sin_metadatos'
cache_key = pa.signature({'dem': pa.file_hash(DEM_PATH), 'codigo': implementacion,
    'versiones': versiones_cache, 'python': sys.version, 'metodo': 'FilledElevation_defaults_D8_Area',
    'aggradation_slope': getattr(d.FilledElevation, 'aggradation_slope', None)})
recuperar = False
if REUTILIZAR_HIDROLOGIA and pa.output('topografia', 'manifest', 'json').exists():
    try:
        previo = pa.read_stage('topografia')
        recuperar = previo.get('cache_key') == cache_key
    except (ValueError, OSError, KeyError) as exc:
        print('Caché no reutilizable:', exc)
pa.begin('topografia')
print('Recuperar hidrología:', recuperar)


# Notebook cell 11
if recuperar:
    filled = d.FilledElevation.load(previo['files']['filled']['path'])
    flow_direction = d.FlowDirectionD8.load(previo['files']['flow_direction']['path'])
    area = d.Area.load(previo['files']['area']['path'])
else:
    elevation = d.Elevation.load(str(DEM_PATH))
    if elevation._griddata.shape != shape_dem:
        raise ValueError('Las dimensiones cargadas no coinciden con el DEM.')

    indices_topo = elevation._xy_to_rowscols(catalogo[['x', 'y']].itertuples(index=False, name=None))
    indices_qc = list(qc_outlets[['fila', 'columna']].itertuples(index=False, name=None))
    if list(indices_topo) != indices_qc:
        raise ValueError('La conversión de coordenadas difiere del control previo; revisar TopoAnalysis.')
    print('DEM original cargado.')
    filled = d.FilledElevation(elevation=elevation)
    print('Relleno terminado.')
    flow_direction = d.FlowDirectionD8(flooded_dem=filled)
    print('Dirección D8 terminada.')
    area = d.Area(flow_direction=flow_direction)
    print('Área aportante terminada.')
    log_area = None
    if GUARDAR_LOGAREA:
        log_area = d.LogArea(area=area)
        print('Logaritmo del área terminado.')
log_area = d.LogArea(area=area) if GUARDAR_LOGAREA else None


# Notebook cell 13
productos = {
    'filled': (filled, DATA_DIR / f'dem_filled_{resolucion_m:g}m_{FECHA}.tif'),
    'flow_direction': (flow_direction, DATA_DIR / f'flow_direction_D8_{resolucion_m:g}m_{FECHA}.tif'),
    'area': (area, DATA_DIR / f'area_aportante_{resolucion_m:g}m_{FECHA}.tif'),
}
if GUARDAR_LOGAREA:
    productos['log_area'] = (log_area, DATA_DIR / f'logarea_{resolucion_m:g}m_{FECHA}.tif')

rutas_productos = {}
for nombre, (grilla, ruta) in productos.items():
    grilla.save(str(ruta))
    with rasterio.open(ruta) as src:
        if (src.crs != crs_dem or src.shape != shape_dem
                or not src.transform.almost_equals(transform_dem)):
            raise ValueError(f'Georreferenciación incompatible: {ruta.name}')
    rutas_productos[nombre] = str(ruta)
    print('Guardado:', ruta.name)

versiones = {}
for paquete in ('TopoAnalysis', 'numpy', 'rasterio', 'geopandas', 'pandas', 'shapely'):
    try:
        versiones[paquete] = metadata.version(paquete)
    except metadata.PackageNotFoundError:
        versiones[paquete] = 'sin metadatos de distribución'

manifiesto = {
    'fecha_ejecucion': FECHA,
    'dem_original': str(DEM_PATH),
    'dem_size_bytes': DEM_PATH.stat().st_size,
    'dem_mtime_ns': DEM_PATH.stat().st_mtime_ns,
    'crs': crs_dem.to_string(),
    'transformacion': list(transform_dem)[:6],
    'dimensiones': list(shape_dem),
    'resolucion_m': resolucion_m,
    'area_celda_m2': area_celda_m2,
    'celdas_validas_dem': n_validas,
    'metodo': 'FilledElevation -> FlowDirectionD8 -> Area',
    'parametros': 'Valores por defecto de la implementación instalada',
    'topoanalysis_fuente': str(Path(d.__file__).resolve()),
    'python': sys.version,
    'ejecutable': sys.executable,
    'versiones': versiones,
    'productos': rutas_productos,
    'catalogo_outlets': str(DATA_DIR / f'catalogo_outlets_{FECHA}.csv'),
    'estado_cuencas': 'pendiente',
}
ruta_manifiesto = DATA_DIR / f'manifiesto_topografico_{FECHA}.json'
ruta_manifiesto.write_text(json.dumps(manifiesto, ensure_ascii=False, indent=2), encoding='utf-8')


# Notebook cell 15
estados = []
ruta_resumen = DATA_DIR / f'resumen_cuencas_{FECHA}.csv'
manifiesto['estado_cuencas'] = 'en_proceso'
ruta_manifiesto.write_text(json.dumps(manifiesto, ensure_ascii=False, indent=2), encoding='utf-8')

for rio in catalogo.itertuples(index=False):
    estado = {'rio_id': rio.rio_id, 'nombre': rio.nombre, 'x': rio.x, 'y': rio.y,
              'estado': 'error', 'area_km2': np.nan, 'n_celdas': 0,
              'mascara': '', 'poligono': '', 'detalle': ''}
    basin_mask = image = positivas = gdf = geometrias = None
    try:
        fila, columna = flow_direction._xy_to_rowscols([(rio.x, rio.y)])[0]
        if fila is None or columna is None:
            raise ValueError('Outlet fuera de la grilla.')
        prefijo = f'rio_{int(rio.rio_id):02d}'
        ruta_mask = MASK_DIR / f'{prefijo}_mask_{FECHA}.tif'
        ruta_poly = POLY_DIR / f'{prefijo}_cuenca_{FECHA}.gpkg'

        basin_mask = d.Mask(flow_direction=flow_direction, outlets=((rio.x, rio.y),))
        basin_mask.save(str(ruta_mask))
        estado['mascara'] = str(ruta_mask)
        basin_mask = None

        with rasterio.open(ruta_mask) as src:
            if src.crs != crs_dem:
                raise ValueError('El CRS de la máscara no coincide con el DEM.')
            if src.shape != shape_dem or not src.transform.almost_equals(transform_dem):
                raise ValueError('La máscara no está alineada con el DEM.')
            image = src.read(1, masked=True)
            positivas = (~np.ma.getmaskarray(image)
                         & np.isfinite(image.data) & (image.data > 0))
            n_celdas = int(np.count_nonzero(positivas))
            if n_celdas == 0:
                raise ValueError('La cuenca está vacía.')
            if not positivas[fila, columna]:
                raise ValueError('La máscara no incluye el outlet.')
            area_m2 = n_celdas * area_celda_m2
            geometrias = [shape(geom) for geom, _ in shapes(
                positivas.astype('uint8'), mask=positivas, transform=src.transform)]

        gdf = gpd.GeoDataFrame({'rio_id': [rio.rio_id] * len(geometrias)},
                               geometry=geometrias, crs=crs_dem)
        if gdf.empty or not gdf.geometry.is_valid.all():
            raise ValueError('Se generaron geometrías vacías o inválidas.')
        if not np.isclose(float(gdf.geometry.area.sum()), area_m2, rtol=1e-8, atol=0.01):
            raise ValueError('El área vectorial no coincide con la máscara.')
        gdf.to_file(ruta_poly, layer='cuenca', driver='GPKG', mode='w')
        estado.update(estado='ok', area_km2=area_m2 / 1e6,
                      n_celdas=n_celdas, poligono=str(ruta_poly))
        print(f'{int(rio.rio_id):02d} | {rio.nombre}: {area_m2 / 1e6:.2f} km²')
    except Exception as exc:
        estado['detalle'] = f'{type(exc).__name__}: {exc}'
        print(f'{int(rio.rio_id):02d} | {rio.nombre}: ERROR — {exc}')
    finally:
        estados.append(estado)
        pd.DataFrame(estados).to_csv(ruta_resumen, index=False, encoding='utf-8-sig')
        basin_mask = image = positivas = gdf = geometrias = None

resumen = pd.DataFrame(estados)
fallidos = resumen.loc[resumen['estado'] != 'ok']
manifiesto.update(
    estado_cuencas='completo' if fallidos.empty else 'incompleto',
    resumen_cuencas=str(ruta_resumen),
    cuencas_correctas=int(resumen['estado'].eq('ok').sum()),
    cuencas_esperadas=len(catalogo),
)
ruta_manifiesto.write_text(json.dumps(manifiesto, ensure_ascii=False, indent=2), encoding='utf-8')
display(resumen[['rio_id', 'nombre', 'estado', 'area_km2', 'detalle']])
print('Resumen:', ruta_resumen)
if not fallidos.empty:
    raise RuntimeError('Hay cuencas pendientes. Consultar la columna detalle del resumen.')


# Notebook cell 16
archivos = {**rutas_productos, 'outlets': DATA_DIR / f'catalogo_outlets_{FECHA}.csv',
            'qc_outlets': DATA_DIR / f'qc_outlets_{FECHA}.csv', 'resumen': ruta_resumen,
            'procedencia': ruta_manifiesto}
for row in resumen.itertuples(index=False):
    archivos[f'mascara_{row.rio_id}'] = row.mascara
    archivos[f'poligono_{row.rio_id}'] = row.poligono
pa.finish('topografia', archivos, extra={'cache_key': cache_key, 'crs': CFG['crs'],
    'dem_original': str(DEM_PATH), 'dem_sha256': pa.file_hash(DEM_PATH)})
print('Topografía completa para esta corrida:', RUN_ID)


# Notebook cell 17

