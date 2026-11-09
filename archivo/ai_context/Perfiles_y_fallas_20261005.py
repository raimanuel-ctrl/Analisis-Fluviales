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

pa.NOTEBOOK_PATH = CARPETA_NOTEBOOKS / 'Perfiles_y_fallas_20261005.ipynb'


# Notebook cell 4
import geopandas as gpd
from shapely.geometry import Point, LineString
from matplotlib.colors import is_color_like
perfiles, nombres, manifest = pa.load_profiles()
manuales, manuales_rechazados = pa.load_manuals(perfiles)
deps = pa.dependency('perfiles')
deps[CFG['kp_manual_csv']] = pa.file_hash(CFG['kp_manual_csv'])
pa.begin('fallas')
geo_crs = CFG['crs']


# Notebook cell 5
def geo_texto(valor):
    return "" if pd.isna(valor) else str(valor).strip()

def geo_exigir_columnas(tabla, columnas, nombre):
    faltantes = set(columnas) - set(tabla.columns)
    if faltantes:
        raise ValueError(f"{nombre}: faltan columnas {sorted(faltantes)}")

def geo_leer_vector(ruta, tipos):
    capa = gpd.read_file(ruta)

    if capa.empty or capa.crs is None:
        raise ValueError(f"Capa vacía o sin CRS definido:\n{ruta}")

    if (
        capa.geometry.isna().any()
        or capa.geometry.is_empty.any()
        or not capa.geometry.is_valid.all()
        or not capa.geom_type.isin(tipos).all()
    ):
        raise ValueError(f"Revisar geometrías de:\n{ruta}")

    # Reproyección en memoria; no modifica el shapefile original.
    return capa.to_crs(geo_crs).reset_index(drop=True)

def geo_partes(geometria):
    """Extraer puntos o líneas de una intersección."""
    if geometria.is_empty:
        return
    if hasattr(geometria, "geoms"):
        for parte in geometria.geoms:
            yield from geo_partes(parte)
    else:
        yield geometria

def geo_posiciones(linea, interseccion):
    """Distancias sobre el cauce de puntos y extremos de intervalos."""
    posiciones = []

    for parte in geo_partes(interseccion):
        if parte.geom_type == "Point":
            posiciones.append(linea.project(parte))

        elif parte.geom_type in ("LineString", "LinearRing"):
            posiciones.extend([
                linea.project(Point(parte.coords[0])),
                linea.project(Point(parte.coords[-1])),
            ])

    return posiciones

def geo_unidad_pendiente(codigo, estado):
    return (
        codigo in ("SIN_ASIG", "Hf/Hl1")
        or "pendiente" in str(estado).lower()
        or estado == "No documentado"
    )


# Notebook cell 6
poligonos = geo_leer_vector(Path(CFG['geologia_integrada']), ['Polygon','MultiPolygon'])
cartas = geo_leer_vector(Path(CFG['fallas_integradas']), ['LineString','MultiLineString'])
chaf = geo_leer_vector(Path(CFG['chaf']), ['LineString','MultiLineString'])
geo_exigir_columnas(poligonos, ['codigo','fuente','estado'], 'Unidades integradas')
geo_exigir_columnas(cartas, ['uid','fuente','falla'], 'Trazas integradas')
geo_exigir_columnas(chaf, ['f_id','falla'], 'CHAF')
for col in ['codigo','fuente','estado','color_hex']:
    if col not in poligonos: poligonos[col] = ''
    poligonos[col] = poligonos[col].fillna('').astype(str)
poligonos['codigo'] = poligonos.codigo.replace('', 'SIN_ASIG')
cartas['id_traza'] = cartas.fuente.astype(str) + ':' + cartas.uid.astype(str)
cartas['uso_metrico'] = False  # Precisión métrica de digitalizaciones no validada.
cartas['nombre'] = cartas.falla.fillna('')
chaf['fuente'], chaf['uso_metrico'] = 'CHAF', True
chaf['nombre'] = chaf.falla.fillna('')
chaf['id_traza'] = 'CHAF:' + chaf.f_id.astype(str)
capas = [cartas, chaf]
fuentes_usadas = ['geologia_integrada','fallas_integradas','chaf']
if CFG['incluir_maule_en_calculo']:
    maule = geo_leer_vector(Path(CFG['maule']), ['LineString','MultiLineString'])
    maule['fuente'], maule['uso_metrico'] = 'Maule: contexto', False
    maule['id_traza'] = ['Maule:'+str(i) for i in range(len(maule))]
    if 'nombre' not in maule: maule['nombre'] = 'Maule'
    capas.append(maule); fuentes_usadas.append('maule')
fallas_integradas = gpd.GeoDataFrame(pd.concat(capas, ignore_index=True), geometry='geometry', crs=geo_crs)
for key in fuentes_usadas:
    path = Path(CFG[key])
    for component in path.parent.glob(path.stem + '.*'):
        if component.is_file(): deps[str(component)] = pa.file_hash(component)


# Notebook cell 8
filas_tramos, filas_contactos, filas_cruces, filas_kp, filas_estado = [], [], [], [], []
resultados_geo = {}
for rid in pa.EXPECTED:
    estado = {'ID': rid, 'rio': nombres[rid], 'estado': 'error', 'detalle': '',
              'cobertura_pct': np.nan, 'registros_cruce': 0, 'KP_validos': int(manuales.ID.eq(rid).sum())}
    try:
        p = perfiles[rid]
        linea = LineString(p[['x_utm','y_utm']].to_numpy(float))
        if not linea.is_simple or linea.length <= 0: raise ValueError('Cauce no simple o longitud nula')
        s, chi = p.distancia_m.to_numpy(float), p.chi_m.to_numpy(float)
        if not np.isclose(linea.length, s[-1], atol=.01, rtol=1e-8): raise ValueError('Longitud geométrica incompatible con perfil')
        cortes = [0., linea.length]
        for i in poligonos.sindex.query(linea, predicate='intersects'):
            cortes += geo_posiciones(linea, linea.intersection(poligonos.geometry.iloc[i].boundary))
        cortes = np.unique(np.clip(np.round(cortes,3),0,linea.length)); cortes[0],cortes[-1] = 0.,linea.length
        tramos = []
        for a,b in zip(cortes[:-1],cortes[1:]):
            if b<=a: continue
            ids = poligonos.sindex.query(linea.interpolate((a+b)/2), predicate='intersects')
            codigo, situacion, fuente = 'SIN_COBERTURA','Sin cobertura',''
            if len(ids)>1: codigo,situacion = 'AMBIGUO','Solape o límite compartido'
            elif len(ids)==1:
                poly = poligonos.iloc[ids[0]]; codigo,situacion,fuente = poly.codigo,poly.estado,poly.fuente
            tramos.append({'ID':rid,'desde_m':a,'hasta_m':b,'codigo':codigo,'estado':situacion,'fuente':fuente})
        contactos = []
        for a,b in zip(tramos[:-1],tramos[1:]):
            if a['codigo']=='SIN_COBERTURA' or b['codigo']=='SIN_COBERTURA': tipo='borde_cobertura'
            elif a['codigo']=='AMBIGUO' or b['codigo']=='AMBIGUO': tipo='limite_ambiguo'
            elif a['fuente'] != b['fuente']: tipo='costura_cartas'
            elif a['codigo'] != b['codigo']: tipo='cambio_unidad_pendiente_revision'
            else: continue
            contactos.append({'ID':rid,'distancia_m':a['hasta_m'],'tipo':tipo,
                'codigo_antes':a['codigo'],'codigo_despues':b['codigo'],'fuente_antes':a['fuente'],'fuente_despues':b['fuente']})
        cruces = []
        for i in fallas_integradas.sindex.query(linea,predicate='intersects'):
            f = fallas_integradas.iloc[i]
            for parte in geo_partes(linea.intersection(f.geometry)):
                pos = geo_posiciones(linea,parte)
                if not pos: continue
                a,b = min(pos),max(pos)
                cruces.append({'ID':rid,'id_traza':f.id_traza,'fragmento':int(i),'fuente':f.fuente,
                    'nombre':f.nombre,'tipo':'puntual' if parte.geom_type=='Point' else 'coincidencia_lineal',
                    'desde_m':a,'hasta_m':b,'chi_desde_m':float(np.interp(a,s,chi)),
                    'chi_hasta_m':float(np.interp(b,s,chi)),'uso_metrico':bool(f.uso_metrico)})
        for k in manuales.loc[manuales.ID.eq(rid)].itertuples():
            punto = Point(k.x_utm,k.y_utm)
            metricas = fallas_integradas.loc[fallas_integradas.uso_metrico]
            row = {'ID':rid,'id_kp':k.id_kp,'perfil_hash':pa.profile_hash(p),
                   'estado':'sin_trazas_metricas','id_traza':'','fuente':'','distancia_plana_m':np.nan,
                   'distancia_al_cruce_sobre_cauce_m':np.nan}
            if not metricas.empty:
                distances = metricas.geometry.distance(punto)
                i = distances.idxmin(); f = metricas.loc[i]
                row.update(estado='proximidad_geometrica',id_traza=f.id_traza,fuente=f.fuente,
                    distancia_plana_m=float(distances.loc[i]))
            candidatos = [max(c['desde_m']-k.distancia_m, k.distancia_m-c['hasta_m'],0.) for c in cruces if c['uso_metrico']]
            if candidatos: row['distancia_al_cruce_sobre_cauce_m']=min(candidatos)
            filas_kp.append(row)
        cobertura = 100*(1-sum(t['hasta_m']-t['desde_m'] for t in tramos if t['codigo']=='SIN_COBERTURA')/linea.length)
        filas_tramos += tramos; filas_contactos += contactos; filas_cruces += cruces
        resultados_geo[rid] = {'tramos':tramos,'contactos':contactos,'cruces':cruces}
        estado.update(estado='ok',cobertura_pct=cobertura,registros_cruce=len(cruces))
    except (ValueError, KeyError, RuntimeError) as exc:
        estado['detalle'] = str(exc)
    filas_estado.append(estado)
estados_geo = pd.DataFrame(filas_estado)
tablas_geo = {
    'estados': estados_geo,
    'tramos': pd.DataFrame(filas_tramos, columns=['ID','desde_m','hasta_m','codigo','estado','fuente']),
    'contactos': pd.DataFrame(filas_contactos, columns=['ID','distancia_m','tipo','codigo_antes','codigo_despues','fuente_antes','fuente_despues']),
    'cruces': pd.DataFrame(filas_cruces, columns=['ID','id_traza','fragmento','fuente','nombre','tipo','desde_m','hasta_m','chi_desde_m','chi_hasta_m','uso_metrico']),
    'kp_fallas': pd.DataFrame(filas_kp, columns=['ID','id_kp','perfil_hash','estado','id_traza','fuente','distancia_plana_m','distancia_al_cruce_sobre_cauce_m']),
    'manuales_rechazados': manuales_rechazados}
archivos = {k:pa.save_table('fallas',k,t) for k,t in tablas_geo.items()}
display(estados_geo)
if not estados_geo.estado.eq('ok').all(): raise RuntimeError('Análisis geológico incompleto; consultar estados.')
pa.finish('fallas', archivos, deps, extra={'fuentes': fuentes_usadas,
    'advertencia':'Fragmentos e intersecciones no equivalen a fallas independientes; cartografía preliminar'})


# Notebook cell 10
RIOS_VER = list(perfiles)
MOSTRAR_FUENTES = None  # None muestra todas; filtra exclusivamente la figura.
GUARDAR_FIGURAS = False
MOSTRAR_TERRAZAS = False
terrazas_vista = pa.stage_table('terrazas','asociaciones') if MOSTRAR_TERRAZAS else pd.DataFrame()
codigos = sorted(set(t['codigo'] for r in resultados_geo.values() for t in r['tramos']))
colores = {c:plt.get_cmap('tab20')(i%20) for i,c in enumerate(codigos)}
colores.update(SIN_COBERTURA='#dddddd', AMBIGUO='#555555', SIN_ASIG='#aaaaaa')
for rid in map(str,RIOS_VER):
    p, r = perfiles[rid], resultados_geo[rid]
    fig,(ax,banda) = plt.subplots(2,1,sharex=True,figsize=(12,6),gridspec_kw={'height_ratios':[10,1]})
    ax.plot(p.chi_m,p.elevacion_m,color='steelblue',lw=1)
    for t in r['tramos']:
        a,b = np.interp([t['desde_m'],t['hasta_m']],p.distancia_m,p.chi_m)
        banda.axvspan(a,b,color=colores[t['codigo']],label=t['codigo'])
    for c in r['contactos']:
        x = np.interp(c['distancia_m'],p.distancia_m,p.chi_m)
        if c['tipo']=='costura_cartas': banda.axvline(x,color='black',ls='--')
        elif c['tipo']=='cambio_unidad_pendiente_revision': ax.axvline(x,color='brown',ls=':',alpha=.5)
    for c in r['cruces']:
        if MOSTRAR_FUENTES is not None and c['fuente'] not in MOSTRAR_FUENTES: continue
        if c['tipo']=='puntual': ax.axvline(c['chi_desde_m'],color='gray',ls='--',alpha=.5)
        else: ax.axvspan(c['chi_desde_m'],c['chi_hasta_m'],color='gray',alpha=.2)
    for k in manuales.loc[manuales.ID.eq(rid)].itertuples():
        x = p.chi_m.iloc[int(k.indice)]
        ax.scatter(x,k.elevacion_m,marker='*',s=100,color='red')
        ax.annotate(k.id_kp,(x,k.elevacion_m),fontsize=7)
    if MOSTRAR_TERRAZAS:
        for t in terrazas_vista.loc[terrazas_vista.ID.eq(rid)].to_dict('records'):
            ax.scatter(t['chi_proyectado_m'],t['Elevación (m)'],marker='^',color='green')
    handles,labels=banda.get_legend_handles_labels(); unique=dict(zip(labels,handles))
    ax.legend(unique.values(),unique.keys(),fontsize=7,bbox_to_anchor=(1.01,1),loc='upper left')
    ax.set(title=nombres[rid]+' — cartografía preliminar',ylabel='Elevación (m)')
    banda.set(xlabel='χ (m)',yticks=[]); fig.tight_layout()
    if GUARDAR_FIGURAS: fig.savefig(pa.output('fallas',f'fig_rio_{int(rid):02d}','png'),dpi=300)
    plt.show();plt.close(fig)


# Notebook cell 12
EJECUTAR_TRANSECTAS = False
MOSTRAR_TRANSECTAS = False


# Notebook cell 13
if EJECUTAR_TRANSECTAS:
    from pathlib import Path
    from datetime import datetime
    import json
    import warnings

    import numpy as np
    import pandas as pd
    import geopandas as gpd
    import rasterio
    import matplotlib.pyplot as plt

    from shapely.geometry import LineString
    from IPython.display import display

    BASE = Path(CFG['base_dir'])
    RUTA_DEM, RUTA_FALLAS = Path(CFG['dem']), Path(CFG['chaf'])
    RUTA_SALIDA = pa.ROOT / 'transectas'
    pa.begin('transectas')
    # ── Parámetros exploratorios ──────────────────────────────────────────
    LONGITUD_SECTOR_M = 2000.0
    DISTANCIA_POR_LADO_M = 1000.0
    SEPARACION_TRANSECTAS_M = 100.0
    VENTANA_RUMBO_LOCAL_M = 200.0

    # Número mínimo de perfiles con elevación válida para resumir
    # una posición transversal.
    MIN_PERFILES_VALIDOS = 3

    # Umbral geométrico: evita asignar bloques cuando la perpendicular
    # local es poco compatible con el sentido de buzamiento registrado.
    MIN_ALINEACION_BUZAMIENTO = 0.5

    for ruta in (RUTA_DEM, RUTA_FALLAS):
        if not ruta.exists():
            raise FileNotFoundError(ruta)

    for valor in (
        LONGITUD_SECTOR_M,
        DISTANCIA_POR_LADO_M,
        SEPARACION_TRANSECTAS_M,
        VENTANA_RUMBO_LOCAL_M,
    ):
        if not np.isfinite(valor) or valor <= 0:
            raise ValueError("Las distancias deben ser positivas y finitas.")

    plt.rcParams.update({
        "figure.figsize": (10, 5),
        "font.size": 10,
        "savefig.dpi": 200,
    })

    print("Configuración preparada.")
    with rasterio.open(RUTA_DEM) as src:
        CRS_DEM = src.crs
        RESOLUCION_M = max(abs(src.res[0]), abs(src.res[1]))
        LIMITES_DEM = src.bounds

    if CRS_DEM is None or CRS_DEM.to_epsg() != 32718:
        raise ValueError("Se esperaba el DEM en EPSG:32718.")

    fallas = gpd.read_file(RUTA_FALLAS)

    if fallas.crs is None:
        raise ValueError("La capa de fallas no tiene CRS.")

    fallas = fallas.to_crs(CRS_DEM)

    requeridas = {"f_id", "falla", "rumbo", "buz", "dir_buz", "sentido"}
    faltantes = requeridas - set(fallas.columns)
    if faltantes:
        raise ValueError(f"Faltan atributos: {faltantes}")

    if (
        fallas.geometry.isna().any()
        or fallas.geometry.is_empty.any()
        or (~fallas.geometry.is_valid).any()
    ):
        raise ValueError("Hay geometrías ausentes, vacías o inválidas.")

    if not fallas.geom_type.isin(["LineString", "MultiLineString"]).all():
        raise ValueError("La capa debe contener únicamente líneas.")

    # No unir registros distintos de una misma falla.
    fallas = fallas.explode(index_parts=False).reset_index(drop=True)
    fallas["registro"] = np.arange(1, len(fallas) + 1)
    fallas["f_id"] = fallas["f_id"].astype(str)

    # Muestreo transversal cercano a la resolución del DEM.
    # Linspace garantiza incluir exactamente el cero.
    n_pasos = int(np.ceil(DISTANCIA_POR_LADO_M / RESOLUCION_M))
    DISTANCIAS = np.linspace(
        -DISTANCIA_POR_LADO_M,
        DISTANCIA_POR_LADO_M,
        2 * n_pasos + 1,
    )

    print(f"Partes de trazas: {len(fallas)}")
    print(f"Resolución del DEM: {RESOLUCION_M:g} m")
    print(f"Ángulos de buzamiento informados: {fallas['buz'].notna().sum()}")
    print("Distancia negativa: lado del buzamiento / bloque colgante.")
    print("Distancia positiva: lado opuesto / bloque yaciente.")

    display(fallas[
        ["registro", "f_id", "falla", "rumbo", "buz", "dir_buz", "sentido"]
    ])


# Notebook cell 14
if EJECUTAR_TRANSECTAS:
    AZIMUT_BUZAMIENTO = {
        "N": 0, "NE": 45, "E": 90, "SE": 135,
        "S": 180, "SW": 225, "W": 270, "NW": 315,
    }

    resultados = {}
    registros_transectas = []
    registros_sectores = []
    incidencias = []

    from functools import lru_cache
    from rasterio.windows import Window
    with rasterio.open(RUTA_DEM) as dem:
        @lru_cache(maxsize=128)
        def leer_bloque(columna, fila, ancho, alto):
            return dem.read(1, window=Window(columna,fila,ancho,alto), masked=True)

        for _, falla in fallas.iterrows():

            linea = falla.geometry
            longitud = linea.length
            registro = int(falla["registro"])
            fid = falla["f_id"]

            direccion = str(falla["dir_buz"]).strip().upper()

            if direccion not in AZIMUT_BUZAMIENTO:
                incidencias.append({
                    "f_id": fid,
                    "sector": None,
                    "motivo": "Sentido de buzamiento no reconocido",
                })
                continue

            if longitud <= 0:
                incidencias.append({
                    "f_id": fid,
                    "sector": None,
                    "motivo": "Longitud nula",
                })
                continue

            azimut = np.radians(AZIMUT_BUZAMIENTO[direccion])
            vector_buz = np.array([np.sin(azimut), np.cos(azimut)])

            # Dividir en sectores similares, evitando un último sector diminuto.
            nsectores = max(1, int(np.ceil(longitud / LONGITUD_SECTOR_M)))
            limites = np.linspace(0, longitud, nsectores + 1)

            for numero in range(nsectores):
                inicio, fin = limites[numero:numero + 2]

                clave = f"{fid}_R{registro:02d}_S{numero + 1:02d}"

                # Muestreo en puntos medios: evita repetir límites de sectores.
                cantidad = max(
                    1,
                    int(np.ceil((fin - inicio) / SEPARACION_TRANSECTAS_M))
                )
                posiciones = (
                    inicio
                    + (np.arange(cantidad) + 0.5) * (fin - inicio) / cantidad
                )

                perfiles_z = []

                for numero_t, posicion in enumerate(posiciones, start=1):
                    mitad = VENTANA_RUMBO_LOCAL_M / 2

                    p0 = linea.interpolate(max(0, posicion - mitad))
                    p1 = linea.interpolate(min(longitud, posicion + mitad))
                    centro = linea.interpolate(posicion)

                    tangente = np.array([p1.x - p0.x, p1.y - p0.y])
                    norma = np.linalg.norm(tangente)

                    if norma == 0:
                        incidencias.append({
                            "f_id": fid, "sector": clave,
                            "motivo": "No se pudo definir la orientación local",
                        })
                        continue

                    tangente /= norma
                    normal = np.array([-tangente[1], tangente[0]])

                    # Normal positiva hacia el buzamiento.
                    if np.dot(normal, vector_buz) < 0:
                        normal *= -1

                    alineacion = float(np.dot(normal, vector_buz))

                    if alineacion < MIN_ALINEACION_BUZAMIENTO:
                        incidencias.append({
                            "f_id": fid, "sector": clave,
                            "motivo": "Orientación local ambigua respecto al buzamiento",
                        })
                        continue

                    # Distancias negativas hacia el buzamiento:
                    # izquierda del gráfico = bloque colgante.
                    xy = (
                        np.array([centro.x, centro.y])[None, :]
                        - DISTANCIAS[:, None] * normal[None, :]
                    )

                    # Muestreo de celda del DEM, sin interpolar ni rellenar vacíos.
                     # Lectura directa por ventanas, sin dem.sample.
                    # Conserva el valor de celda, sin interpolar.
                    from rasterio.windows import Window

                    z = np.full(len(xy), np.nan, dtype=float)

                    # Coordenadas UTM -> posiciones de píxel.
                    inversa = ~dem.transform
                    columnas_float = (
                        inversa.a * xy[:, 0]
                        + inversa.b * xy[:, 1]
                        + inversa.c
                    )
                    filas_float = (
                        inversa.d * xy[:, 0]
                        + inversa.e * xy[:, 1]
                        + inversa.f
                    )

                    dentro = (
                        np.isfinite(columnas_float)
                        & np.isfinite(filas_float)
                        & (columnas_float >= 0)
                        & (columnas_float < dem.width)
                        & (filas_float >= 0)
                        & (filas_float < dem.height)
                    )

                    indices = np.flatnonzero(dentro)

                    if len(indices):
                        columnas = np.floor(
                            columnas_float[indices]
                        ).astype(np.int64)

                        filas = np.floor(
                            filas_float[indices]
                        ).astype(np.int64)

                        # Agrupar puntos en ventanas pequeñas de 256 × 256.
                        # Evita cargar todo el DEM o una ventana muy grande.
                        tamano = 256
                        bloques = np.column_stack([
                            filas // tamano,
                            columnas // tamano,
                        ])

                        bloques_unicos, grupo = np.unique(
                            bloques, axis=0, return_inverse=True
                        )

                        for numero_bloque, (bloque_f, bloque_c) in enumerate(
                            bloques_unicos
                        ):
                            seleccion = grupo == numero_bloque

                            fila_inicio = int(bloque_f * tamano)
                            columna_inicio = int(bloque_c * tamano)

                            alto = min(tamano, dem.height - fila_inicio)
                            ancho = min(tamano, dem.width - columna_inicio)

                            ventana = leer_bloque(columna_inicio, fila_inicio, ancho, alto)

                            valores = ventana[
                                filas[seleccion] - fila_inicio,
                                columnas[seleccion] - columna_inicio,
                            ]

                            z[indices[seleccion]] = (
                                np.ma.asarray(valores, dtype=float)
                                .filled(np.nan)
                            )

                    z[~np.isfinite(z)] = np.nan
                    perfiles_z.append(z)

                    rumbo_local = float(
                        np.degrees(np.arctan2(tangente[0], tangente[1])) % 180
                    )
                    rumbo_atributo = pd.to_numeric(
                        falla["rumbo"], errors="coerce"
                    )
                    diferencia = (
                        abs((rumbo_local - rumbo_atributo + 90) % 180 - 90)
                        if pd.notna(rumbo_atributo) else np.nan
                    )

                    registros_transectas.append({
                        "sector": clave,
                        "transecta": numero_t,
                        "f_id": fid,
                        "falla": falla["falla"],
                        "dir_buz": direccion,
                    "buz": pd.to_numeric(falla.get("buz"), errors="coerce"),
                        "cinematica": falla["sentido"],
                        "rumbo_loc": rumbo_local,
                        "dif_rumbo": diferencia,
                        "alineacion": alineacion,
                        "fr_valid": float(np.isfinite(z).mean()),
                        "geometry": LineString(xy),
                    })

                if not perfiles_z:
                    incidencias.append({
                        "f_id": fid, "sector": clave,
                        "motivo": "Sector sin transectas utilizables",
                    })
                    continue

                matriz = np.vstack(perfiles_z)
                nvalidos = np.isfinite(matriz).sum(axis=0)

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    p10, mediana, p90 = np.nanpercentile(
                        matriz, [10, 50, 90], axis=0
                    )

                suficiente = nvalidos >= MIN_PERFILES_VALIDOS
                p10[~suficiente] = np.nan
                mediana[~suficiente] = np.nan
                p90[~suficiente] = np.nan

                tabla = pd.DataFrame({
                    "distancia_m": DISTANCIAS,
                    "p10_m": p10,
                    "mediana_m": mediana,
                    "p90_m": p90,
                    "n_validos": nvalidos,
                })

                resultados[clave] = {
                    "f_id": fid,
                    "falla": falla["falla"],
                    "tipo": falla.get("tipo", ""),
                    "cinematica": falla["sentido"],
                    "dir_buz": direccion,
                    "buz": pd.to_numeric(falla.get("buz"), errors="coerce"),
                    "z": matriz,
                    "tabla": tabla,
                }

                registros_sectores.append({
                    "sector": clave,
                    "f_id": fid,
                    "falla": falla["falla"],
                    "tipo": falla.get("tipo", ""),
                    "cinematica": falla["sentido"],
                    "inicio_m": inicio,
                    "fin_m": fin,
                    "n_perfiles": len(matriz),
                    "fraccion_valida": float(np.isfinite(matriz).mean()),
                    "fraccion_resumible": float(suficiente.mean()),
                })

            print(f"Revisada: {fid} — {falla['falla']}")

    pa.save_table('transectas', 'incidencias', pd.DataFrame(incidencias, columns=['f_id','sector','motivo']))
    if not resultados:
        raise ValueError('No se generaron transectas; consultar incidencias exportadas.')

    transectas = gpd.GeoDataFrame(
        registros_transectas, geometry="geometry", crs=CRS_DEM
    )
    resumen = pd.DataFrame(registros_sectores)
    incidencias_df = pd.DataFrame(incidencias)

    print(f"\nSectores procesados: {len(resultados)}")
    print(f"Transectas generadas: {len(transectas)}")
    print(f"Incidencias: {len(incidencias_df)}")

    display(resumen)

    if not incidencias_df.empty:
        display(incidencias_df)
    leer_bloque.cache_clear()
    archivos_transectas = {
        'resumen': pa.save_table('transectas','resumen',resumen),
        'incidencias': pa.save_table('transectas','incidencias',incidencias_df)}
    ruta_transectas = pa.output('transectas','transectas','gpkg')
    transectas.to_file(ruta_transectas,layer='transectas',driver='GPKG',mode='w')
    archivos_transectas['geometrias'] = ruta_transectas
    for clave,r in resultados.items():
        safe = re.sub(r'[^A-Za-z0-9_]+','_',str(clave))
        archivos_transectas['sector_'+safe] = pa.save_table('transectas','sector_'+safe,r['tabla'])
        archivos_transectas['individuales_'+safe] = pa.save_table('transectas','individuales_'+safe,
            pd.DataFrame(r['z'],columns=[f'{x:.6f}m' for x in DISTANCIAS]))
    deps_transectas = {str(RUTA_DEM):pa.file_hash(RUTA_DEM)}
    for component in RUTA_FALLAS.parent.glob(RUTA_FALLAS.stem+'.*'):
        if component.is_file(): deps_transectas[str(component)] = pa.file_hash(component)
    pa.finish('transectas',archivos_transectas,deps_transectas,extra={'parametros':{
        'longitud_sector_m':LONGITUD_SECTOR_M,'distancia_por_lado_m':DISTANCIA_POR_LADO_M,
        'separacion_transectas_m':SEPARACION_TRANSECTAS_M,'ventana_rumbo_m':VENTANA_RUMBO_LOCAL_M,
        'min_perfiles_validos':MIN_PERFILES_VALIDOS,'alineacion_minima':MIN_ALINEACION_BUZAMIENTO}})


# Notebook cell 15
if EJECUTAR_TRANSECTAS and MOSTRAR_TRANSECTAS:
    # ── VISUALIZACIÓN: solo cortes, con orientación geográfica ────────────

    SECTORES_VER = (
        resumen.groupby("f_id", sort=False)["sector"].first().tolist()
    )

    # Para elegir sectores concretos:
    SECTORES_VER = ["300904_R03_S01", "300903_R02_S01"]

    MOSTRAR_INDIVIDUALES = True
    MOSTRAR_BANDA = False


    def nombre_direccion(vector):
        """Dirección aproximada de un vector Este–Norte."""
        azimut = np.degrees(np.arctan2(vector[0], vector[1])) % 360

        nombres = [
            "N · Norte", "NE · Noreste",
            "E · Este", "SE · Sureste",
            "S · Sur", "SW · Suroeste",
            "W · Oeste", "NW · Noroeste",
        ]
        return nombres[int(np.floor((azimut + 22.5) / 45)) % 8]


    def figura_sector(clave):
        r = resultados[clave]
        tabla = r["tabla"]
        t = transectas[transectas["sector"] == clave]

        if t.empty:
            raise ValueError(f"No hay transectas para {clave}.")

        # Vector de cada transecta desde su extremo negativo al positivo.
        vectores = np.array([
            np.asarray(g.coords[-1]) - np.asarray(g.coords[0])
            for g in t.geometry
        ], dtype=float)

        vectores /= np.linalg.norm(vectores, axis=1)[:, None]
        vector_medio = vectores.mean(axis=0)

        norma = np.linalg.norm(vector_medio)
        if norma < 0.5:
            raise ValueError(
                f"{clave}: orientaciones demasiado variables. "
                "Conviene revisar o reducir el sector."
            )

        vector_medio /= norma

        # Queremos el extremo más al norte a la izquierda.
        # Para cortes prácticamente E–W, colocamos el oeste a la izquierda.
        invertir = (
            vector_medio[1] > 0
            if abs(vector_medio[1]) > 1e-6
            else vector_medio[0] < 0
        )

        signo = -1 if invertir else 1

        x_original = tabla["distancia_m"].to_numpy()
        x_grafico = signo * x_original
        orden = np.argsort(x_grafico)
        x = x_grafico[orden]

        z = r["z"][:, orden]
        mediana = tabla["mediana_m"].to_numpy()[orden]
        p10 = tabla["p10_m"].to_numpy()[orden]
        p90 = tabla["p90_m"].to_numpy()[orden]
        n_validos = tabla["n_validos"].to_numpy()[orden]

        vector_derecha = signo * vector_medio
        izquierda = nombre_direccion(-vector_derecha)
        derecha = nombre_direccion(vector_derecha)

        bloque_izq = "Lado opuesto al buzamiento" if invertir else "Lado del buzamiento"
        bloque_der = "Lado del buzamiento" if invertir else "Lado opuesto al buzamiento"

        fig, ax = plt.subplots(figsize=(11, 5.5))

        if MOSTRAR_INDIVIDUALES:
            ax.plot(
                x, z.T,
                color="gray", linewidth=0.7, alpha=0.3
            )
            ax.plot(
                [], [], color="gray", linewidth=0.8,
                label="Perfiles individuales"
            )

        if MOSTRAR_BANDA:
            ax.fill_between(
                x, p10, p90,
                color="steelblue", alpha=0.2,
                label="80 % central de las elevaciones"
            )

        ax.plot(
            x, mediana,
            color="navy", linewidth=2.2,
            label="Mediana de los perfiles"
        )

        # ── Plano de falla esquemático ───────────────────────────────────
        buz = pd.to_numeric(r.get('buz'), errors='coerce')
        MANTEO_SUPUESTO = float(buz) if pd.notna(buz) and 0 < buz <= 90 else 60.0
        origen_buz = 'informado en CHAF' if pd.notna(buz) and 0 < buz <= 90 else 'supuesto'
        ax.set_aspect('equal', adjustable='datalim')
        PROFUNDIDAD_DIBUJADA_M = 100.0  # Solo controla la extensión visual

        # Anclar en la mediana del relieve, exactamente en la traza.
        indice_cero = int(np.argmin(np.abs(x)))
        z_traza = mediana[indice_cero]

        if np.isfinite(z_traza):
            # En la convención original, el colgante tiene distancia negativa.
            # "signo" incorpora la inversión geográfica de la figura.
            sentido_horizontal = -signo

            profundidad = np.linspace(0, PROFUNDIDAD_DIBUJADA_M, 50)
            desplazamiento = (
                profundidad / np.tan(np.radians(MANTEO_SUPUESTO))
            )

            x_plano = sentido_horizontal * desplazamiento
            z_plano = z_traza - profundidad

            ax.plot(
                x_plano, z_plano,
                color="black",
                linestyle="--",
                linewidth=2,
                label=f"Plano esquemático: {MANTEO_SUPUESTO:g}° ({origen_buz})",
            )

            ax.scatter(
                [0], [z_traza],
                color="black", s=35, zorder=5,
                label="Traza en superficie",
            )

        else:
            # No inventar una elevación para anclar el plano.
            ax.text(
                0.5, 0.03,
                "Sin mediana válida en la traza: no se dibuja el plano.",
                transform=ax.transAxes,
                ha="center", fontsize=9, color="darkred",
            )

        # Marcar posiciones donde faltan perfiles válidos para la mediana.
        insuficiente = n_validos < MIN_PERFILES_VALIDOS
        ax.fill_between(
            x, 0, 1,
            where=insuficiente,
            transform=ax.get_xaxis_transform(),
            color="darkorange", alpha=0.12,
            label="Datos insuficientes para la mediana"
            if insuficiente.any() else None,
        )

        ax.text(
            0.01, 1.02,
            f"{izquierda}\n{bloque_izq}",
            transform=ax.transAxes,
            ha="left", va="bottom", fontweight="bold",
        )
        ax.text(
            0.99, 1.02,
            f"{derecha}\n{bloque_der}",
            transform=ax.transAxes,
            ha="right", va="bottom", fontweight="bold",
        )

        ax.set(
            xlim=(x.min(), x.max()),
            xlabel=(
                "Distancia transversal a la falla (m)\n"
                "0 = traza cartografiada; valores negativos a la izquierda"
            ),
            ylabel="Elevación del DEM (m)",
        )

        fig.suptitle(
            f"{r['falla']} — {clave}\n"
            f"{r['cinematica']} | {r['tipo']} | "
            f"Buzamiento hacia {r['dir_buz']}",
            y=1.04,
        )

        ax.grid(alpha=0.2)
        ax.legend(loc="best", fontsize=9)
        fig.tight_layout()

        return fig


    for clave in SECTORES_VER:
        if clave not in resultados:
            raise KeyError(f"Sector inexistente: {clave}")

        fig = figura_sector(clave)
        plt.show()
        plt.close(fig)


# Notebook cell 16

