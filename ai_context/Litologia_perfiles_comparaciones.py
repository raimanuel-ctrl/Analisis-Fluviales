# Code-only reading snapshot; not an independently executable script.
# Notebook is authoritative. Cell magics, external data and cell order still apply.

# Notebook cell 2
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import os
import sys
import json
import hashlib

# Ayuda a localizar las DLL de Conda en Windows.
# Debe ejecutarse antes de importar NumPy/SciPy.
DLL_HANDLES = []
carpeta_dll = Path(sys.prefix) / "Library" / "bin"

if os.name == "nt" and carpeta_dll.exists():
    os.environ["PATH"] = (
        str(carpeta_dll) + os.pathsep + os.environ.get("PATH", "")
    )
    DLL_HANDLES.append(os.add_dll_directory(str(carpeta_dll)))

for variable, subcarpeta in [
    ("GDAL_DATA", "gdal"),
    ("PROJ_DATA", "proj"),
]:
    carpeta = Path(sys.prefix) / "Library" / "share" / subcarpeta
    if carpeta.exists():
        os.environ[variable] = str(carpeta)

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.stats import linregress
from scipy.spatial import cKDTree
from IPython.display import display

# ---------------- MODIFICAR AQUÍ ----------------

ETAPA = "litologia"   
BASE = Path(
    r"C:\datos"
    r"\Memoria\Codigos Stanford"
)

RUN = "20261005_post_auditoria_01"
PRODUCTOS = BASE / "Productos" / "post_auditoria" / RUN

THETA = 0.40
A0 = 1e6                 # m²
N = 1.0

# None carga los 18 ríos.
# Para cargar solamente los seleccionados:
# IDS_CARGAR = ["3", "4", "7", "8", "9", "10"]
IDS_CARGAR = None

RIOS_K = ["3", "4", "7", "8", "9", "10"]

GUARDAR_TABLAS = True
GUARDAR_FIGURAS = False
GUARDAR_PDF = True

FECHA = datetime.now(
    ZoneInfo("America/Los_Angeles")
).strftime("%Y%m%d_%H%M%S_%f")

SALIDA = BASE / "Productos" / ETAPA / FECHA


FUENTES = {}


def registrar_archivo(ruta):
    """Registrar la huella de un archivo leído."""
    ruta = Path(ruta)
    huella = hashlib.sha256(ruta.read_bytes()).hexdigest()
    FUENTES[str(ruta.resolve())] = huella
    return huella


def leer_manifest(etapa):
    ruta = PRODUCTOS / etapa / f"manifest_{RUN}.json"
    registrar_archivo(ruta)

    documento = json.loads(ruta.read_text(encoding="utf-8-sig"))

    if documento.get("estado") != "completo":
        raise ValueError(f"La etapa {etapa} no está completa.")

    if documento.get("run_id") != RUN:
        raise ValueError(f"La etapa {etapa} pertenece a otra corrida.")

    return documento


def leer_csv_manifest(manifest, clave):
    item = manifest["files"][clave]
    ruta = Path(item["path"])

    if registrar_archivo(ruta) != item["sha256"]:
        raise ValueError(
            f"El archivo cambió después de exportarse: {ruta.name}"
        )

    return pd.read_csv(
        ruta,
        dtype={"ID": str},
        float_precision="round_trip",
    )


MANIFEST_PERFILES = leer_manifest("perfiles")

if MANIFEST_PERFILES["config"].get("crs") != "EPSG:32718":
    raise ValueError("Verificar el CRS antes de continuar.")

ids_disponibles = sorted(MANIFEST_PERFILES["rios"], key=int)

IDS = (
    ids_disponibles
    if IDS_CARGAR is None
    else list(map(str, IDS_CARGAR))
)

if set(IDS) - set(ids_disponibles):
    raise ValueError("Hay IDs que no existen en la exportación.")

NOMBRES = {
    rid: MANIFEST_PERFILES["rios"][rid]["nombre"]
    for rid in IDS
}

perfiles = {}
filas_qc = []

for rid in IDS:
    original = leer_csv_manifest(
        MANIFEST_PERFILES, f"rio_{rid}"
    )

    columnas = [
        "distancia_m",
        "elevacion_m",
        "area_m2",
        "x_utm",
        "y_utm",
    ]

    if set(columnas) - set(original):
        raise ValueError(f"Río {rid}: faltan columnas.")

    # Conservar únicamente geometría y datos originales necesarios.
    p = original[columnas].copy()

    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    a = p.area_m2.to_numpy(float)

    if (
        len(p) < 2
        or not np.isfinite(s).all()
        or not np.isclose(s[0], 0)
        or np.any(np.diff(s) <= 0)
    ):
        raise ValueError(f"Río {rid}: distancia inválida.")

    if not np.isfinite(a).all() or np.any(a <= 0):
        raise ValueError(f"Río {rid}: área inválida.")

    if not np.isfinite(
        p[["x_utm", "y_utm"]].to_numpy(float)
    ).all():
        raise ValueError(f"Río {rid}: coordenadas inválidas.")

    if not np.isfinite(z[0]):
        raise ValueError(f"Río {rid}: elevación del outlet inválida.")

    p["chi_040_m"] = cumulative_trapezoid(
        (A0 / a) ** THETA,
        x=s,
        initial=0,
    )

    perfiles[rid] = p

    filas_qc.append({
        "ID": rid,
        "rio": NOMBRES[rid],
        "nodos": len(p),
        "longitud_km": s[-1] / 1000,
        "z_outlet_m": z[0],
        "huecos_z": int((~np.isfinite(z)).sum()),
    })

qc_perfiles = pd.DataFrame(filas_qc)

display(qc_perfiles)
print("Carpeta de resultados:", SALIDA)

# Notebook cell 4
RUTA_MANUALES = BASE / "kp_manuales_28_09.csv"

INCLUIR_CANDIDATOS = True
TOLERANCIA_XY_M = 0.10
TOLERANCIA_Z_M = 0.10

aceptados = []
rechazados = []

# ---------------- Selecciones manuales ----------------

if RUTA_MANUALES.exists():
    registrar_archivo(RUTA_MANUALES)

    manuales_originales = pd.read_csv(
        RUTA_MANUALES,
        sep=";",
        decimal=",",
        dtype={"ID": str},
    )

    for _, r in manuales_originales.iterrows():
        rid = str(r["ID"])

        if rid not in perfiles:
            continue

        try:
            p = perfiles[rid]
            xy = np.array(
                [float(r["x_utm"]), float(r["y_utm"])]
            )

            if not np.isfinite(xy).all():
                raise ValueError("Coordenadas no finitas.")

            distancia_xy, indice = cKDTree(
                p[["x_utm", "y_utm"]].to_numpy()
            ).query(xy)

            indice = int(indice)
            nodo = p.iloc[indice]

            diferencia_z = abs(
                float(nodo.elevacion_m) - float(r["elevacion_m"])
            )

            if (
                distancia_xy > TOLERANCIA_XY_M
                or diferencia_z > TOLERANCIA_Z_M
            ):
                raise ValueError(
                    f"Desajuste XY={distancia_xy:.3f} m; "
                    f"Z={diferencia_z:.3f} m."
                )

            aceptados.append({
                "ID": rid,
                "id_kp": str(r["id_kp"]),
                "tipo": "manual_reubicado",
                "indice": indice,
                "distancia_m": float(nodo.distancia_m),
                "elevacion_m": float(nodo.elevacion_m),
                "chi_040_m": float(nodo.chi_040_m),
                "x_utm": float(nodo.x_utm),
                "y_utm": float(nodo.y_utm),
                "theta_origen": r.get("theta_seleccion", np.nan),
                "indice_original": r["indice"],
                "distancia_original_m": r["distancia_m"],
                "hash_original": r["perfil_hash"],
                "hash_actual": MANIFEST_PERFILES["rios"][rid][
                    "perfil_hash"
                ],
                "desajuste_xy_m": distancia_xy,
                "desajuste_z_m": diferencia_z,
            })

        except (ValueError, TypeError, KeyError) as exc:
            rechazados.append({
                "ID": rid,
                "id_kp": r.get("id_kp", ""),
                "origen": "manual",
                "motivo": str(exc),
            })
else:
    print("No se encontró el CSV de selecciones manuales.")

# ---------------- Candidatos exportados ----------------

ruta_manifest_analisis = (
    PRODUCTOS
    / "analisis_perfiles"
    / f"manifest_{RUN}.json"
)

if INCLUIR_CANDIDATOS and ruta_manifest_analisis.exists():
    manifest_analisis = leer_manifest("analisis_perfiles")

    candidatos = leer_csv_manifest(
        manifest_analisis, "catalogo_candidatos"
    )

    for _, r in candidatos.iterrows():
        rid = str(r["ID"])

        if rid not in perfiles:
            continue

        # Los manuales se incorporan desde su archivo de selección.
        if str(r.get("metodo", "")).lower() == "manual":
            continue

        try:
            p = perfiles[rid]
            valor_indice = float(r["indice"])

            if (
                not np.isfinite(valor_indice)
                or not valor_indice.is_integer()
                or not 0 <= valor_indice < len(p)
            ):
                raise ValueError("Índice inválido.")

            indice = int(valor_indice)
            nodo = p.iloc[indice]

            if r["perfil_hash"] != MANIFEST_PERFILES["rios"][rid][
                "perfil_hash"
            ]:
                raise ValueError("Huella de perfil incompatible.")

            for columna in [
                "distancia_m",
                "elevacion_m",
                "x_utm",
                "y_utm",
            ]:
                if not np.isclose(
                    float(r[columna]),
                    float(nodo[columna]),
                    rtol=0,
                    atol=0.01,
                ):
                    raise ValueError(f"Desajuste en {columna}.")

            aceptados.append({
                "ID": rid,
                "id_kp": str(r["id_kp"]),
                "tipo": "candidato",
                "indice": indice,
                "distancia_m": float(nodo.distancia_m),
                "elevacion_m": float(nodo.elevacion_m),
                "chi_040_m": float(nodo.chi_040_m),
                "x_utm": float(nodo.x_utm),
                "y_utm": float(nodo.y_utm),
                "theta_origen": r.get("theta", np.nan),
            })

        except (ValueError, TypeError, KeyError) as exc:
            rechazados.append({
                "ID": rid,
                "id_kp": r.get("id_kp", ""),
                "origen": "candidato",
                "motivo": str(exc),
            })

columnas_kp = [
    "ID", "id_kp", "tipo", "indice",
    "distancia_m", "elevacion_m", "chi_040_m",
    "x_utm", "y_utm", "theta_origen",
    "indice_original", "distancia_original_m",
    "hash_original", "hash_actual",
    "desajuste_xy_m", "desajuste_z_m",
]

kp = pd.DataFrame(aceptados, columns=columnas_kp)

kp = kp.drop_duplicates(
    ["ID", "id_kp", "tipo"]
).reset_index(drop=True)

# Si un candidato coincide exactamente con un nodo manual,
# conservar la selección manual y evitar duplicar el símbolo.
claves_manuales = set(
    zip(
        kp.loc[kp.tipo.eq("manual_reubicado"), "ID"],
        kp.loc[kp.tipo.eq("manual_reubicado"), "indice"],
    )
)

duplicado_manual = pd.Series(
    [
        (rid, indice) in claves_manuales
        for rid, indice in zip(kp.ID, kp.indice)
    ],
    index=kp.index,
    dtype=bool,
)

kp = kp.loc[
    ~(kp.tipo.eq("candidato") & duplicado_manual)
].reset_index(drop=True)

kp_rechazados = pd.DataFrame(
    rechazados,
    columns=["ID", "id_kp", "origen", "motivo"],
)

display(kp)
display(kp_rechazados)

print(
    "Manuales reubicados:",
    int(kp.tipo.eq("manual_reubicado").sum()),
)
print("Candidatos:", int(kp.tipo.eq("candidato").sum()))

# Notebook cell 6
# LITOLOGÍA — INTERSECCIÓN Y MISMA UNIDAD EN AMBAS RIBERAS
import geopandas as gpd

GEO_DIR = Path(
    r"C:\datos\Memoria_Codex\outputs"
    r"\Geologia_integrada_Chanco_Jorquera_Gajardo"
)

CRS = "EPSG:32718"

for nombre in ["Geologia_unidades", "Geologia_fallas"]:
    for extension in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
        archivo = GEO_DIR / f"{nombre}{extension}"
        if archivo.exists():
            registrar_archivo(archivo)

geo = gpd.read_file(GEO_DIR / "Geologia_unidades.shp")
fallas = gpd.read_file(GEO_DIR / "Geologia_fallas.shp")

if geo.crs is None or fallas.crs is None:
    raise ValueError("Alguna capa no tiene CRS definido.")

geo = geo.to_crs(CRS)
fallas = fallas.to_crs(CRS)

for nombre, capa in [("geología", geo), ("fallas", fallas)]:
    if (
        capa.geometry.isna().any()
        or capa.geometry.is_empty.any()
        or not capa.geometry.is_valid.all()
    ):
        raise ValueError(f"Revisar geometrías de {nombre}.")

CAMPOS_GEO = [
    "codigo", "unidad", "grupo_lito", "fuente",
    "color_hex", "uso_metr", "escala",
]

if set(CAMPOS_GEO) - set(geo.columns):
    raise ValueError("Faltan campos esperados en la geología.")


# REEMPLAZAR SOLO LA FUNCIÓN asignar_unidades DEL BLOQUE DE LITOLOGÍA.
# geo y CAMPOS_GEO deben estar cargados antes de esta definición.
from shapely.geometry import Point, LineString

# Límite operativo para vacíos estrechos; editable, no precisión de la carta.
MAX_DIST_RIBERA_M = 20.0
VENTANA_DIRECCION_M = 30.0
RIOS_RELLENO_RIBERAS = {str(i) for i in range(7, 19)}
RIOS_SIN_CARTA_CONFIRMADA = {"3", "4"}

def _unidad_en_ribera(punto, normal, signo):
    extremo = Point(punto.x + signo*normal[0]*MAX_DIST_RIBERA_M,
                    punto.y + signo*normal[1]*MAX_DIST_RIBERA_M)
    rayo = LineString([punto, extremo])
    indices = list(geo.sindex.intersection(rayo.bounds))
    candidatas = geo.iloc[indices].copy()
    candidatas = candidatas.loc[candidatas.geometry.intersects(rayo)].copy()
    if candidatas.empty:
        return None, np.nan
    intersecciones = candidatas.geometry.intersection(rayo)
    distancias = intersecciones.distance(punto)
    minimo = float(distancias.min())
    cercanas = candidatas.loc[np.isclose(distancias, minimo, rtol=0, atol=1e-6)]
    if len(cercanas[["codigo", "fuente"]].drop_duplicates()) != 1:
        return None, minimo     # empate entre unidades: no forzar asignación
    return cercanas.iloc[0], minimo

def asignar_unidades(puntos, ids, posiciones_m):
    """Intersección directa; vacíos solo si ambas riberas tienen la misma unidad/fuente."""
    if MAX_DIST_RIBERA_M <= 0 or VENTANA_DIRECCION_M <= 0:
        raise ValueError("Las distancias de búsqueda deben ser positivas.")
    ids = [str(ids)]*len(puntos) if isinstance(ids, (str, int)) else list(map(str, ids))
    posiciones = np.asarray(posiciones_m, dtype=float)
    if len(ids) != len(puntos) or len(posiciones) != len(puntos):
        raise ValueError("IDs, posiciones y puntos tienen longitudes diferentes.")
    if not puntos.index.is_unique:
        raise ValueError("Los índices de puntos deben ser únicos.")
    union = gpd.sjoin(puntos[["geometry"]], geo[CAMPOS_GEO+["geometry"]],
                      how="left", predicate="intersects")
    grupos = {indice: grupo for indice, grupo in union.groupby(level=0, sort=False)}
    filas = []
    for indice, rid, posicion, punto in zip(puntos.index, ids, posiciones, puntos.geometry):
        registro = {campo: None for campo in CAMPOS_GEO}
        registro.update(indice=indice, estado_geo="Sin cobertura", etiqueta="Sin cobertura",
                        color_hex="#B7B7B7", metodo_asignacion="Sin asignar",
                        distancia_poligono_m=np.nan, distancia_izquierda_m=np.nan,
                        distancia_derecha_m=np.nan, motivo_asignacion="")
        if rid in RIOS_SIN_CARTA_CONFIRMADA:
            registro["motivo_asignacion"] = "Sin carta disponible para este análisis"
            filas.append(registro); continue
        validos = grupos[indice].loc[grupos[indice].index_right.notna()]
        if len(validos) == 1:
            r = validos.iloc[0]
            registro.update({campo: r[campo] for campo in CAMPOS_GEO})
            registro.update(estado_geo="asignado", etiqueta=f"{r.codigo} · {r.fuente}",
                            metodo_asignacion="Intersección directa", distancia_poligono_m=0.)
        elif len(validos) > 1:
            registro.update(estado_geo="Ambiguo", etiqueta="Ambiguo",
                            motivo_asignacion="Solape de polígonos: conservar para revisión")
        elif rid not in RIOS_RELLENO_RIBERAS:
            registro["motivo_asignacion"] = "Fuera de los ríos habilitados para rellenar vacíos"
        else:
            p = perfiles[rid]
            s = p.distancia_m.to_numpy(float)
            if not np.isfinite(posicion) or not s[0] <= posicion <= s[-1]:
                raise ValueError(f"Río {rid}: posición fuera del perfil.")
            a, b = max(s[0], posicion-VENTANA_DIRECCION_M), min(s[-1], posicion+VENTANA_DIRECCION_M)
            dx = float(np.interp(b, s, p.x_utm)-np.interp(a, s, p.x_utm))
            dy = float(np.interp(b, s, p.y_utm)-np.interp(a, s, p.y_utm))
            longitud = np.hypot(dx, dy)
            if longitud <= 1e-8:
                registro["motivo_asignacion"] = "Dirección local del cauce indeterminada"
                filas.append(registro); continue
            normal = np.array([-dy, dx])/longitud
            izquierda, di = _unidad_en_ribera(punto, normal, 1)
            derecha, dd = _unidad_en_ribera(punto, normal, -1)
            registro.update(distancia_izquierda_m=di, distancia_derecha_m=dd)
            if izquierda is None or derecha is None:
                registro["motivo_asignacion"] = "No hay dos riberas inequívocas dentro del límite"
            elif (str(izquierda.codigo), str(izquierda.fuente)) != (str(derecha.codigo), str(derecha.fuente)):
                registro.update(estado_geo="Ambiguo", etiqueta="Ambiguo",
                                motivo_asignacion="Unidades o fuentes distintas a ambos lados")
            else:
                r = izquierda
                registro.update({campo: r[campo] for campo in CAMPOS_GEO})
                registro.update(estado_geo="asignado", etiqueta=f"{r.codigo} · {r.fuente}",
                                metodo_asignacion="Misma unidad en ambas riberas",
                                distancia_poligono_m=min(di, dd))
        filas.append(registro)
    return pd.DataFrame(filas).set_index("indice").reindex(puntos.index)



segmentos = {}
cambios = []
propuestas = []

for rid, p in perfiles.items():
    q = pd.DataFrame({
        "s0": p.distancia_m.to_numpy()[:-1],
        "s1": p.distancia_m.to_numpy()[1:],
        "z0": p.elevacion_m.to_numpy()[:-1],
        "z1": p.elevacion_m.to_numpy()[1:],
        "chi0": p.chi_040_m.to_numpy()[:-1],
        "chi1": p.chi_040_m.to_numpy()[1:],
    })

    puntos = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(
            (
                p.x_utm.to_numpy()[:-1]
                + p.x_utm.to_numpy()[1:]
            ) / 2,
            (
                p.y_utm.to_numpy()[:-1]
                + p.y_utm.to_numpy()[1:]
            ) / 2,
        ),
        crs=CRS,
    )

    q = pd.concat([
        q,
        asignar_unidades(puntos, ids=rid,
                         posiciones_m=(q.s0.to_numpy()+q.s1.to_numpy())/2),
    ], axis=1)
    q["ID"] = rid

    q["familia"] = q.grupo_lito.where(
        q.estado_geo.eq("asignado"),
        q.estado_geo,
    ).fillna("Sin cobertura")

    segmentos[rid] = q

    for j in range(1, len(q)):
        anterior = q.iloc[j - 1]
        actual = q.iloc[j]

        if anterior.etiqueta == actual.etiqueta:
            continue

        if (
            anterior.estado_geo != "asignado"
            or actual.estado_geo != "asignado"
        ):
            tipo = "cobertura_o_solape"
        elif anterior.fuente != actual.fuente:
            tipo = "cambio_fuente_posible_costura"
        else:
            tipo = "cambio_unidad_muestreado"

        cambios.append({
            "ID": rid,
            "distancia_m": actual.s0,
            "antes": anterior.etiqueta,
            "despues": actual.etiqueta,
            "tipo": tipo,
        })

    aluvial = (
        q.estado_geo.eq("asignado")
        & q.grupo_lito.fillna("").str.contains(
            "fluvial|aluvial",
            case=False,
            regex=True,
        )
    )

    grupos = (
        aluvial.ne(aluvial.shift())
        | q.etiqueta.ne(q.etiqueta.shift())
    ).cumsum()

    for _, grupo in q.loc[aluvial].groupby(grupos[aluvial]):
        propuestas.append({
            "ID": rid,
            "inicio_m": float(grupo.s0.iloc[0]),
            "fin_m": float(grupo.s1.iloc[-1]),
            "unidad": grupo.etiqueta.iloc[0],
            "grupo": grupo.grupo_lito.iloc[0],
            "estado": "propuesta_cartografica",
        })

cambios_geo = pd.DataFrame(
    cambios,
    columns=["ID", "distancia_m", "antes", "despues", "tipo"],
)

propuestas_aluviales = pd.DataFrame(
    propuestas,
    columns=[
        "ID", "inicio_m", "fin_m",
        "unidad", "grupo", "estado",
    ],
)

# Añadir litología en la ubicación de cada KP.
kp_litologia = kp.copy()

if not kp.empty:
    puntos_kp = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(kp.x_utm, kp.y_utm),
        crs=CRS,
    )

    puntos_kp.index = kp.index
    kp_litologia = pd.concat([
        kp,
        asignar_unidades(puntos_kp, ids=kp.ID,
                         posiciones_m=kp.distancia_m.to_numpy()),
    ], axis=1)

display(propuestas_aluviales)
display(kp_litologia)

# Notebook cell 8
VENTANAS_M = [250, 500, 1000]
THETAS_EXPLORAR = np.arange(0.25, 0.651, 0.025)
MIN_NODOS_PENDIENTE = 5


def pendientes_ventana(s, z, ventana):
    lo = np.searchsorted(s, s - ventana / 2)
    hi = np.searchsorted(
        s, s + ventana / 2, side="right"
    )

    validos = np.isfinite(z)
    zz = np.where(validos, z, 0.0)

    def suma_intervalos(valores):
        acumulada = np.r_[0.0, np.cumsum(valores)]
        return acumulada[hi] - acumulada[lo]

    n = hi - lo
    sx = suma_intervalos(s)
    sy = suma_intervalos(zz)

    varianza_x = suma_intervalos(s * s) - sx * sx / n
    covarianza = suma_intervalos(s * zz) - sx * sy / n

    ok = (
        (suma_intervalos(validos.astype(float)) == n)
        & (n >= MIN_NODOS_PENDIENTE)
        & (varianza_x > 0)
    )

    pendiente = np.full(len(s), np.nan)
    pendiente[ok] = covarianza[ok] / varianza_x[ok]

    # Control independiente en hasta 100 ventanas.
    indices_control = np.unique(
        np.linspace(
            0, len(s) - 1, min(100, len(s))
        ).astype(int)
    )

    for j in indices_control:
        if not ok[j]:
            continue

        referencia = linregress(
            s[lo[j]:hi[j]] - s[j],
            z[lo[j]:hi[j]],
        ).slope

        if not np.isclose(
            pendiente[j],
            referencia,
            rtol=1e-5,
            atol=1e-8,
        ):
            raise ValueError(
                f"Control de pendiente falló: "
                f"nodo {j}, ventana {ventana} m."
            )

    return pendiente


sensibilidad_ksn = {}
sensibilidad_theta = {}
filas_resumen = []

for rid, p in perfiles.items():
    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    a = p.area_m2.to_numpy(float)

    sensibilidad_ksn[rid] = {}

    for ventana in VENTANAS_M:
        pendiente = pendientes_ventana(s, z, ventana)

        ksn = np.where(
            np.isfinite(pendiente) & (pendiente > 0),
            pendiente * a**THETA,
            np.nan,
        )

        sensibilidad_ksn[rid][ventana] = ksn

        if ventana == 500:
            p["pendiente_500"] = pendiente
            p["ksn_040"] = ksn

    filas_theta = []

    for theta in THETAS_EXPLORAR:
        chi = cumulative_trapezoid(
            (A0 / a)**theta,
            x=s,
            initial=0,
        )

        r2 = (
            linregress(chi, z).rvalue**2
            if np.isfinite(z).all()
            else np.nan
        )

        filas_theta.append({
            "theta": float(theta),
            "R2": float(r2),
        })

    sensibilidad_theta[rid] = pd.DataFrame(filas_theta)

    r2_referencia = (
        linregress(p.chi_040_m, z).rvalue**2
        if np.isfinite(z).all()
        else np.nan
    )

    ksn_validos = p.ksn_040.dropna()

    filas_resumen.append({
        "ID": rid,
        "rio": NOMBRES[rid],
        "longitud_km": s[-1] / 1000,
        "R2_chi_040": r2_referencia,
        "ksn_mediana_040": (
            float(ksn_validos.median())
            if len(ksn_validos)
            else np.nan
        ),
        "ksn_invalidos": int(p.ksn_040.isna().sum()),
        "manuales": int(
            ((kp.ID == rid) & kp.tipo.eq("manual_reubicado")).sum()
        ),
        "candidatos": int(
            ((kp.ID == rid) & kp.tipo.eq("candidato")).sum()
        ),
    })

resumen_rios = pd.DataFrame(filas_resumen)
display(resumen_rios)

# Notebook cell 9

#OPCIONES GRÁFICAS


get_ipython().run_line_magic("matplotlib", "inline")

import matplotlib
import matplotlib.pyplot as plt

from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages

# ---------------- MODIFICAR AQUÍ ----------------

RIOS_VER = list(perfiles)
# RIOS_VER = ["3", "4", "7", "8", "9", "10"]
# RIOS_VER = ["9"]

COLOR_POR = "familia"       # "familia" o "unidad"
GROSOR = 2.0
BORDE_NEGRO = True
RELLENO = False
BANDA = False
VERTICALES = False
ETIQUETAS_KP = False
GRILLA = True

# GROSOR = 3.0
# BORDE_NEGRO = False
# RELLENO = True
# BANDA = True
# VERTICALES = True
# ETIQUETAS_KP = True
# COLOR_POR = "unidad"

# ------------------------------------------------

plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
})

familias = sorted(
    set(pd.concat(segmentos.values()).familia)
)

paleta_familias = {
    nombre: matplotlib.colors.to_hex(
        matplotlib.colormaps["tab20"](i % 20)
    )
    for i, nombre in enumerate(familias)
}

paleta_familias.update({
    "Sin cobertura": "#B7B7B7",
    "Ambiguo": "#444444",
})


def colores_segmentos(q):
    if COLOR_POR == "familia":
        return [
            paleta_familias.get(v, "#B7B7B7")
            for v in q.familia
        ]

    return [
        v if isinstance(v, str)
        and matplotlib.colors.is_color_like(v)
        else "#B7B7B7"
        for v in q.color_hex
    ]


def marcar_kp(ax, rid, eje="distancia", relativo=False):
    p = perfiles[rid]
    referencia = p.elevacion_m.iloc[0] if relativo else 0

    for tipo, marcador, color in [
        ("manual_reubicado", "o", "black"),
        ("candidato", "^", "#A45720"),
    ]:
        q = kp[(kp.ID == rid) & kp.tipo.eq(tipo)]

        x = (
            q.distancia_m / 1000
            if eje == "distancia"
            else q.chi_040_m / 1000
        )

        ax.scatter(
            x,
            q.elevacion_m - referencia,
            marker=marcador,
            facecolors="white",
            edgecolors=color,
            s=24,
            zorder=5,
        )

        if ETIQUETAS_KP:
            for xx, (_, r) in zip(x, q.iterrows()):
                ax.annotate(
                    str(int(r.indice)),
                    (xx, r.elevacion_m - referencia),
                    xytext=(3, 5),
                    textcoords="offset points",
                    fontsize=6,
                )


def dibujar_perfil(ax, rid, eje="distancia", relativo=False):
    p = perfiles[rid]
    q = segmentos[rid]

    referencia = (
        float(p.elevacion_m.iloc[0]) if relativo else 0
    )

    x0 = (
        q.s0.to_numpy() / 1000
        if eje == "distancia"
        else q.chi0.to_numpy() / 1000
    )

    x1 = (
        q.s1.to_numpy() / 1000
        if eje == "distancia"
        else q.chi1.to_numpy() / 1000
    )

    z0 = q.z0.to_numpy() - referencia
    z1 = q.z1.to_numpy() - referencia

    lineas = np.stack(
        [np.c_[x0, z0], np.c_[x1, z1]],
        axis=1,
    )

    validos = np.isfinite(lineas).all(axis=(1, 2))
    colores = np.asarray(colores_segmentos(q))

    if BORDE_NEGRO:
        ax.add_collection(LineCollection(
            lineas[validos],
            colors="black",
            linewidths=GROSOR + 1.2,
        ))

    ax.add_collection(LineCollection(
        lineas[validos],
        colors=colores[validos],
        linewidths=GROSOR,
    ))

    base = float(p.elevacion_m.iloc[0]) - referencia

    if RELLENO:
        poligonos = [
            [(a, base), (a, z), (b, w), (b, base)]
            for a, b, z, w in zip(
                x0[validos], x1[validos],
                z0[validos], z1[validos],
            )
        ]

        ax.add_collection(PolyCollection(
            poligonos,
            facecolors=colores[validos],
            edgecolors="none",
            alpha=0.20,
        ))

    if BANDA:
        poligonos = [
            [(a, 0), (a, 0.04), (b, 0.04), (b, 0)]
            for a, b in zip(x0, x1)
        ]

        ax.add_collection(PolyCollection(
            poligonos,
            facecolors=colores,
            edgecolors="none",
            transform=ax.get_xaxis_transform(),
        ))

    if VERTICALES:
        cambios_rio = cambios_geo[
            (cambios_geo.ID == rid)
            & cambios_geo.tipo.eq("cambio_unidad_muestreado")
        ]

        for s in cambios_rio.distancia_m:
            x = (
                s / 1000
                if eje == "distancia"
                else np.interp(
                    s, p.distancia_m, p.chi_040_m
                ) / 1000
            )

            ax.axvline(
                x, color="gray", lw=0.6,
                ls=":", alpha=0.5,
            )

    ax.autoscale_view()
    ax.set_xlim(left=0)

    ax.set(
        xlabel=(
            "Distancia desde outlet (km)"
            if eje == "distancia"
            else "χ (km)"
        ),
        ylabel=(
            "Elevación sobre outlet (m)"
            if relativo
            else "Elevación (m)"
        ),
    )

    marcar_kp(ax, rid, eje, relativo)
    ax.grid(GRILLA, alpha=0.16)


def dibujar_ksn(ax, rid):
    p = perfiles[rid]

    for ventana, color in [
        (250, "#AFC4CD"),
        (500, "#174F68"),
        (1000, "#D87929"),
    ]:
        ax.plot(
            p.distancia_m / 1000,
            sensibilidad_ksn[rid][ventana],
            lw=0.9,
            color=color,
            label=f"{ventana} m",
        )

    ax.set(
        xlabel="Distancia desde outlet (km)",
        ylabel=r"$k_{sn}$ (m$^{0.8}$)",
        xlim=(0, None),
    )
    ax.legend(ncol=3)
    ax.grid(GRILLA, alpha=0.16)


def dibujar_pendiente_area(ax, rid):
    p = perfiles[rid]

    validos = (
        np.isfinite(p.pendiente_500)
        & (p.pendiente_500 > 0)
    )

    ax.scatter(
        p.loc[validos, "area_m2"] / 1e6,
        p.loc[validos, "pendiente_500"],
        s=2,
        alpha=0.35,
        color="#2A687D",
        rasterized=True,
    )

    ax.set(
        xscale="log",
        yscale="log",
        xlabel="Área (km²)",
        ylabel="Pendiente S (m/m)",
    )
    ax.grid(GRILLA, alpha=0.16)


def dibujar_theta(ax, rid):
    q = sensibilidad_theta[rid]

    ax.plot(
        q.theta, q.R2,
        "o-", ms=3, lw=1,
        color="#2A687D",
    )

    ax.axvline(THETA, color="#D87929", ls="--", lw=1)

    ax.set(
        xlabel="Concavidad impuesta θ",
        ylabel="R² del perfil z-χ completo",
    )
    ax.grid(GRILLA, alpha=0.16)


def dibujar_mapa(ax, rid):
    p = perfiles[rid]

    xmin = p.x_utm.min() - 500
    xmax = p.x_utm.max() + 500
    ymin = p.y_utm.min() - 500
    ymax = p.y_utm.max() + 500

    g = geo.cx[xmin:xmax, ymin:ymax]
    f = fallas.cx[xmin:xmax, ymin:ymax]

    if not g.empty:
        g.plot(
            ax=ax,
            color=[
                paleta_familias.get(str(v), "#DDDDDD")
                for v in g.grupo_lito
            ],
            edgecolor="none",
            alpha=0.25,
            rasterized=True,
        )

    if not f.empty:
        f.plot(
            ax=ax,
            color="#A13E56",
            linewidth=0.6,
            alpha=0.65,
        )

    ax.plot(p.x_utm, p.y_utm, color="#153B4A", lw=1.2)

    for tipo, marcador in [
        ("manual_reubicado", "o"),
        ("candidato", "^"),
    ]:
        q = kp[(kp.ID == rid) & kp.tipo.eq(tipo)]

        ax.scatter(
            q.x_utm, q.y_utm,
            marker=marcador,
            facecolors="white",
            edgecolors="black",
            s=20,
            zorder=5,
        )

    ax.scatter(
        [p.x_utm.iloc[0]],
        [p.y_utm.iloc[0]],
        marker="s",
        color="#D87929",
        s=20,
        zorder=6,
    )

    ax.set(
        xlim=(xmin, xmax),
        ylim=(ymin, ymax),
        xlabel="Este UTM (m)",
        ylabel="Norte UTM (m)",
    )

    ax.set_aspect("equal")
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.xaxis.set_major_locator(
        matplotlib.ticker.MaxNLocator(3)
    )
    ax.yaxis.set_major_locator(
        matplotlib.ticker.MaxNLocator(4)
    )


def leyenda_litologia(fig, ids):
    q = pd.concat([segmentos[rid] for rid in ids])

    if COLOR_POR == "familia":
        etiquetas = sorted(set(q.familia))
        colores = [paleta_familias[v] for v in etiquetas]
    else:
        q = q.drop_duplicates("etiqueta")
        etiquetas = q.etiqueta.tolist()
        colores = colores_segmentos(q)

    handles = [
        Line2D([0], [0], color=color, lw=3, label=etiqueta)
        for etiqueta, color in zip(etiquetas, colores)
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=min(4, max(1, len(handles))),
        fontsize=6,
        frameon=False,
    )


def cerrar_figura(fig, nombre, pdf=None):
    if GUARDAR_FIGURAS:
        SALIDA.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            SALIDA / f"{nombre}_{FECHA}.png",
            dpi=250,
            bbox_inches="tight",
        )

    if pdf is not None:
        pdf.savefig(fig, dpi=170)

    plt.show()
 #  plt.close(fig)

# Notebook cell 11
# ============================================================
# L4 — SELECTOR MODULAR DE FIGURAS (CONTROL TOTAL)
# ============================================================
get_ipython().run_line_magic("matplotlib", "inline")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Los seis gráficos de cada ficha.
# Las funciones se definen en la celda anterior: # OPCIONES GRÁFICAS.

PANELES = [
    (
        "Distancia–elevación y litología",
        lambda ax, rid: dibujar_perfil(
            ax, rid, eje="distancia", relativo=False
        ),
    ),
    (
        "χ–elevación y litología",
        lambda ax, rid: dibujar_perfil(
            ax, rid, eje="chi", relativo=False
        ),
    ),
    (
        "ksn: sensibilidad a la ventana",
        dibujar_ksn,
    ),
    (
        "Pendiente–área",
        dibujar_pendiente_area,
    ),
    (
        "Sensibilidad a la concavidad",
        dibujar_theta,
    ),
    (
        "Mapa del río y geología",
        dibujar_mapa,
    ),
]

# ------------------------------------------------------------
# 1. ¿QUÉ RÍOS QUIERES VER HOY?
# ------------------------------------------------------------
# Pon aquí los IDs de los ríos que quieres ver.
# Ejemplos:
#   RIOS_A_VER = ["9", "10"]   -> Muestra solo Curanilahue y Curanipe
#   RIOS_A_VER = ["3", "4"]    -> Muestra solo Honda y Malhueco
#   RIOS_A_VER = []            -> Si queda vacío, muestra TODOS los ríos
RIOS_A_VER = ["9", "10"]

# ------------------------------------------------------------
# 2. ¿QUÉ TIPOS DE GRÁFICOS QUIERES ACTIVAR? (True / False)
# ------------------------------------------------------------
VER_FICHAS_INDIVIDUALES = True       # Ficha de 6 paneles de los ríos seleccionados
VER_ORIGEN_COMUN_VS_REAL = True      # Origen común (x=0) y origen continental (Este UTM)
VER_LAMINAS_COMPARATIVAS = False     # Comparaciones de a 6 ríos por panel (Atlas)
VER_MAPA_REGIONAL = False            # Mapa conjunto con fallas y litología

# Si quieres que además exporte a PDF lo que se genere
EXPORTAR_PDF_ACTUAL = False
# ------------------------------------------------------------

pdf = None
if GUARDAR_PDF and EXPORTAR_PDF_ACTUAL:
    SALIDA.mkdir(parents=True, exist_ok=True)
    pdf = PdfPages(SALIDA / f"laminas_seleccion_{FECHA}.pdf")

rios_activos = [str(r) for r in RIOS_A_VER] if RIOS_A_VER else list(perfiles.keys())

try:
    # A. FICHAS INDIVIDUALES (6 paneles por río elegido)
    if VER_FICHAS_INDIVIDUALES:
        for rid in rios_activos:
            if rid not in perfiles:
                continue
            fig, axes = plt.subplots(2, 3, figsize=(11.69, 8.27))
            fig.subplots_adjust(left=0.07, right=0.96, bottom=0.15, top=0.87, wspace=0.36, hspace=0.45)
            fig.suptitle(f"Río {rid} · {NOMBRES.get(rid, '')} · θ={THETA:.2f}", fontsize=14, fontweight="bold")
            
            for ax, (titulo, funcion) in zip(axes.flat, PANELES):
                funcion(ax, rid)
                ax.set_title(titulo, loc="left", fontsize=9)
            
            leyenda_litologia(fig, [rid])
            plt.show()
            if pdf is not None:
                pdf.savefig(fig, dpi=170)
            plt.close(fig)

    # B. COMPARACIÓN MULTI-RÍO: ORIGEN COMÚN vs ORIGEN REAL CONTINENTAL
    if VER_ORIGEN_COMUN_VS_REAL and len(rios_activos) > 0:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
        
        for rid in rios_activos:
            if rid not in perfiles:
                continue
            p = perfiles[rid]
            nom = NOMBRES.get(rid, f"Río {rid}")
            
            # 1. Origen Común (0 km en el outlet)
            ax1.plot(p.distancia_m / 1000.0, p.elevacion_m, label=f"{rid} · {nom}", lw=1.8)
            
            # 2. Origen Real Continental (coordenada Este UTM en km)
            ax2.plot(p.x_utm / 1000.0, p.elevacion_m, label=f"{rid} · {nom}", lw=1.8)

        ax1.set_title("Origen Común (Distancia desde Outlet)", fontsize=11, fontweight="bold")
        ax1.set_xlabel("Distancia a lo largo del cauce (km)", fontsize=10)
        ax1.set_ylabel("Elevación (m)", fontsize=10)
        ax1.grid(True, linestyle=":", alpha=0.6)
        ax1.legend(fontsize=8, loc="upper left")

        ax2.set_title("Origen Real (Posición Continental hacia el Este)", fontsize=11, fontweight="bold")
        ax2.set_xlabel("Coordenada Este UTM 18S (km)", fontsize=10)
        ax2.grid(True, linestyle=":", alpha=0.6)

        plt.suptitle(f"Comparación de Cauces Seleccionados ({len(rios_activos)} ríos)", fontsize=13, fontweight="bold")
        plt.tight_layout()
        plt.show()
        if pdf is not None:
            pdf.savefig(fig, dpi=170)
        plt.close(fig)

    # C. LÁMINAS COMPARATIVAS POR TIPO DE GRÁFICO (Estilo Atlas completo)
    if VER_LAMINAS_COMPARATIVAS:
        for numero, (titulo, funcion) in enumerate(PANELES, 1):
            for inicio in range(0, len(rios_activos), 6):
                grupo = rios_activos[inicio:inicio + 6]
                fig, axes = plt.subplots(2, 3, figsize=(11.69, 8.27))
                fig.subplots_adjust(left=0.07, right=0.96, bottom=0.15, top=0.87, wspace=0.36, hspace=0.45)
                fig.suptitle(f"Comparación · {titulo}", fontsize=14)

                for ax, rid in zip(axes.flat, grupo):
                    funcion(ax, rid)
                    ax.set_title(f"{rid} · {NOMBRES.get(rid, '')}", loc="left", fontsize=8)

                for ax in list(axes.flat)[len(grupo):]:
                    ax.axis("off")

                if numero in [1, 2]:
                    leyenda_litologia(fig, grupo)

                plt.show()
                if pdf is not None:
                    pdf.savefig(fig, dpi=170)
                plt.close(fig)

    # D. MAPA REGIONAL CONJUNTO (Geología + Fallas + Knickpoints)
    if VER_MAPA_REGIONAL:
        fig, ax = plt.subplots(figsize=(9, 9))
        xy = np.concatenate([perfiles[rid][["x_utm", "y_utm"]].to_numpy() for rid in rios_activos if rid in perfiles])
        xmin, ymin = xy.min(axis=0) - 2000
        xmax, ymax = xy.max(axis=0) + 2000

        g = geo.cx[xmin:xmax, ymin:ymax]
        f = fallas.cx[xmin:xmax, ymin:ymax]

        if not g.empty:
            g.plot(ax=ax, color=[paleta_familias.get(str(v), "#DDDDDD") for v in g.grupo_lito],
                   edgecolor="none", alpha=0.25, rasterized=True)
        if not f.empty:
            f.plot(ax=ax, color="#A13E56", linewidth=0.8, alpha=0.7)

        for rid in rios_activos:
            p = perfiles[rid]
            ax.plot(p.x_utm, p.y_utm, lw=1.5, label=f"{rid} · {NOMBRES.get(rid, '')}")
            ax.annotate(rid, (p.x_utm.iloc[0], p.y_utm.iloc[0]), xytext=(3, 4),
                        textcoords="offset points", fontsize=8, fontweight="bold")

        ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), xlabel="Este UTM 18S (m)", ylabel="Norte UTM 18S (m)",
               title="Mapa de Cauces, Litología y Fallas")
        ax.set_aspect("equal")
        ax.ticklabel_format(style="plain", useOffset=False)
        ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1))
        plt.tight_layout()
        plt.show()
        if pdf is not None:
            pdf.savefig(fig, dpi=170)
        plt.close(fig)

finally:
    if pdf is not None:
        pdf.close()

# Notebook cell 12
# ============================================================
# CUADRÍCULA 2x2: DISTANCIA & CHI (ORIGEN COMÚN vs CONTINENTAL)
# ============================================================
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ------------------------------------------------------------
# OPCIONES DE USUARIO (EDITAR AQUÍ)
# ------------------------------------------------------------
# Ríos que quieres ver (si queda [] o None, grafica los 18 ríos)
RIOS_A_COMPARAR = [] 

# Control de suavizado
SUAVIZAR = True
VENTANA_SUAVIZADO_M = 500    # Ventana en metros (ej: 250, 500, 1000)

# Paleta continua de Norte a Sur ("turbo", "viridis", "coolwarm_r", "plasma")
COLORMAP = "turbo"
# ------------------------------------------------------------

rios_plot = [str(r) for r in RIOS_A_COMPARAR] if RIOS_A_COMPARAR else list(perfiles.keys())

# Detectar columna de Chi (chi_040_m o chi_m)
primer_p = next(iter(perfiles.values()))
chi_col = "chi_040_m" if "chi_040_m" in primer_p.columns else "chi_m"

# Ordenar los ríos de Norte a Sur según el Y del outlet
y_outlets = {rid: perfiles[rid]["y_utm"].iloc[0] for rid in rios_plot if rid in perfiles}
rios_ordenados = sorted(y_outlets.keys(), key=lambda r: y_outlets[r], reverse=True)

y_vals = list(y_outlets.values())
y_min, y_max = min(y_vals), max(y_vals)
cmap = plt.get_cmap(COLORMAP)

# Figura 2x2
fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharey=True)
(ax11, ax12), (ax21, ax22) = axes

line_handles = []
line_labels = []

for rid in rios_ordenados:
    p = perfiles[rid]
    nombre = NOMBRES.get(rid, f"Río {rid}")
    y0 = y_outlets[rid]
    
    # Color según latitud (1.0 = Norte, 0.0 = Sur)
    norm_pos = (y0 - y_min) / (y_max - y_min) if y_max > y_min else 0.5
    color = cmap(norm_pos)
    
    s_km = p["distancia_m"].to_numpy() / 1000.0
    z_raw = p["elevacion_m"].to_numpy()
    chi_km = p[chi_col].to_numpy() / 1000.0
    x_outlet_km = p["x_utm"].iloc[0] / 1000.0
    
    # Suavizado móvil
    if SUAVIZAR:
        win_nodos = max(3, int(VENTANA_SUAVIZADO_M / 12.0))
        z_plot = pd.Series(z_raw).rolling(window=win_nodos, center=True, min_periods=1).mean().to_numpy()
        # Perfil crudo tenue de fondo
        for ax, x_d in [(ax11, s_km), (ax12, x_outlet_km + s_km), (ax21, chi_km), (ax22, x_outlet_km + chi_km)]:
            ax.plot(x_d, z_raw, color=color, alpha=0.15, lw=0.7)
    else:
        z_plot = z_raw
        
    line, = ax11.plot(s_km, z_plot, color=color, lw=1.9)
    ax12.plot(x_outlet_km + s_km, z_plot, color=color, lw=1.9)
    ax21.plot(chi_km, z_plot, color=color, lw=1.9)
    ax22.plot(x_outlet_km + chi_km, z_plot, color=color, lw=1.9)
    
    line_handles.append(line)
    line_labels.append(f"{rid} · {nombre}")
    
    # --------------------------------------------------------
    # SOLO KP MANUALES (Se eliminaron completamente los candidatos)
    # --------------------------------------------------------
    if "kp" in globals() and isinstance(kp, pd.DataFrame):
        # Filtramos estrictamente por tipo manual
        filtro_manual = kp["tipo"].astype(str).str.contains("manual", case=False, na=False)
        q = kp[(kp["ID"].astype(str) == str(rid)) & filtro_manual]
        
        for _, row in q.iterrows():
            s_kp = row["distancia_m"] / 1000.0
            z_kp = row["elevacion_m"]
            chi_kp = (row.get("chi_040_m", row.get("chi_m", 0))) / 1000.0
            
            # Círculo con borde negro marcado en los 4 paneles
            for ax, x_kp in [(ax11, s_kp), (ax12, x_outlet_km + s_kp), 
                             (ax21, chi_kp), (ax22, x_outlet_km + chi_kp)]:
                ax.scatter([x_kp], [z_kp], marker="o", color=color, edgecolors="black",
                           s=45, linewidths=1.2, zorder=6)

# Elemento ficticio para indicar el símbolo de KP en la leyenda
kp_dummy = plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="gray", 
                      markeredgecolor="black", markersize=8, label="KP manual")

# --- Configuración Fila Superior (Distancia) ---
ax11.set_title("A. Distancia vs. Elevación (Origen Común)", fontsize=11, fontweight="bold")
ax11.set_xlabel("Distancia a lo largo del cauce (km)", fontsize=10)
ax11.set_ylabel("Elevación (m)", fontsize=10)
ax11.grid(True, linestyle=":", alpha=0.6)

ax12.set_title("B. Distancia Continental (Este Outlet + s)", fontsize=11, fontweight="bold")
ax12.set_xlabel("Coordenada Este desplazada (km)", fontsize=10)
ax12.grid(True, linestyle=":", alpha=0.6)

# --- Configuración Fila Inferior (Chi) ---
ax21.set_title(f"C. χ vs. Elevación (Origen Común, θ={THETA:.2f})", fontsize=11, fontweight="bold")
ax21.set_xlabel("Coordenada χ (km)", fontsize=10)
ax21.set_ylabel("Elevación (m)", fontsize=10)
ax21.grid(True, linestyle=":", alpha=0.6)

ax22.set_title(f"D. χ Continental (Este Outlet + χ, θ={THETA:.2f})", fontsize=11, fontweight="bold")
ax22.set_xlabel("Coordenada Este desplazada en χ (km)", fontsize=10)
ax22.grid(True, linestyle=":", alpha=0.6)

# --- Leyenda Externa a la Derecha (Sin estorbar las curvas) ---
fig.legend(line_handles + [kp_dummy], line_labels + ["KP manual"],
           loc="center left", bbox_to_anchor=(0.905, 0.5), fontsize=8,
           title="Norte → Sur", title_fontsize=9.5, frameon=True)

suav_txt = f"Suavizado: {VENTANA_SUAVIZADO_M}m" if SUAVIZAR else "Sin suavizar"
plt.suptitle(f"Comparación Fluvial Regional ({len(rios_ordenados)} ríos) · {suav_txt}", fontsize=13, fontweight="bold")
plt.subplots_adjust(left=0.06, right=0.89, top=0.93, bottom=0.07, wspace=0.12, hspace=0.25)
plt.show()



# Notebook cell 13
# CELDA NUEVA AL FINAL DE Litologia_perfiles_comparaciones.ipynb
# Ejecutar después de actualizar segmentos y kp_litologia.
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
from matplotlib.lines import Line2D
from matplotlib.colors import is_color_like
from shapely.geometry import box
import matplotlib.pyplot as plt

RIOS_REGIONALES = []     # Vacío: todos los cargados. Ejemplo: ["3", "4", "7", "8", "9", "10"]
GUARDAR_REGIONAL = False
RUTA_CONTEXTO = Path(
    r"C:\datos\Memoria_Codex\outputs"
    r"\Comparacion_KP_fallas_20261004\Resumen_comparacion.md"
)

# Recuperar las etiquetas regionales ya revisadas; no asignar bloques nuevos.
contexto = {}
if RUTA_CONTEXTO.exists():
    registrar_archivo(RUTA_CONTEXTO)
    for linea in RUTA_CONTEXTO.read_text(encoding="utf-8-sig").splitlines():
        if not linea.startswith("|"):
            continue
        c = [v.strip() for v in linea.split("|")[1:-1]]
        if len(c) == 6 and c[0].isdigit():
            contexto[c[0]] = c[2]

ids_region = sorted(map(str, perfiles), key=int) if not RIOS_REGIONALES else list(map(str, RIOS_REGIONALES))
if set(ids_region)-set(perfiles):
    raise ValueError("Algún río solicitado no está cargado.")
filas = []
for rid in ids_region:
    manuales = kp_litologia.loc[
        kp_litologia.ID.astype(str).eq(rid) & kp_litologia.tipo.eq("manual_reubicado")
    ].drop_duplicates("id_kp")
    q = segmentos[rid]
    longitudes = q.s1.to_numpy(float)-q.s0.to_numpy(float)
    asignado = q.estado_geo.eq("asignado").to_numpy()
    unidades, metodos = [], []
    for _, r in manuales.iterrows():
        unidades.append(f"{r.codigo} · {r.fuente}" if r.estado_geo == "asignado" else str(r.estado_geo))
        metodos.append(str(r.get("metodo_asignacion", r.estado_geo)))
    hay_rechazo = False
    if not kp_rechazados.empty:
        hay_rechazo = (kp_rechazados.ID.astype(str).eq(rid) & kp_rechazados.origen.eq("manual")).any()
    filas.append({
        "ID": rid, "Río": NOMBRES[rid], "Calibración K": "Sí" if rid in RIOS_K else "Regional",
        "KP manual": str(len(manuales)) if len(manuales) else ("Rechazado: revisar" if hay_rechazo else "Sin selección"),
        "Elevación KP (m)": " / ".join(f"{v:.1f}" for v in manuales.elevacion_m) or "—",
        "Unidad en KP": " / ".join(unidades) or "—",
        "Método litológico KP": " / ".join(metodos) or "—",
        "% recorrido asignado": 100*longitudes[asignado].sum()/longitudes.sum(),
        "Posición regional (hipótesis)": contexto.get(rid, "Sin KP seleccionado" if manuales.empty else "Pendiente de revisión"),
    })
tabla_regional = pd.DataFrame(filas)
display(tabla_regional.style.format({"% recorrido asignado":"{:.1f}"}))
print("El porcentaje corresponde a longitud de segmentos asignados, no a lecho observado.")
print("Las etiquetas estructurales provienen de la revisión regional; no demuestran causa del KP.")

# Mapa sencillo: todos los cauces, seis ríos de K, KP manuales y fallas disponibles.
xmin = min(float(perfiles[r].x_utm.min()) for r in ids_region)-3000
xmax = max(float(perfiles[r].x_utm.max()) for r in ids_region)+3000
ymin = min(float(perfiles[r].y_utm.min()) for r in ids_region)-3000
ymax = max(float(perfiles[r].y_utm.max()) for r in ids_region)+3000
ventana = box(xmin, ymin, xmax, ymax)
fig, ax = plt.subplots(figsize=(9, 12), constrained_layout=True)
fondo = geo.loc[geo.geometry.intersects(ventana)].copy()
if not fondo.empty:
    fondo.geometry = fondo.geometry.intersection(ventana)
    colores = [c if isinstance(c, str) and is_color_like(c) else "#cccccc" for c in fondo.color_hex]
    fondo.plot(ax=ax, color=colores, edgecolor="0.6", linewidth=.3, alpha=.3)
trazas = fallas.loc[fallas.geometry.intersects(ventana)].copy()
if not trazas.empty:
    trazas.geometry = trazas.geometry.intersection(ventana)
    trazas.plot(ax=ax, color="#bd392f", linewidth=1.1, zorder=3)
for rid in ids_region:
    p = perfiles[rid]
    ax.plot(p.x_utm, p.y_utm, color="#146fa1" if rid in RIOS_K else "0.4",
            lw=1.8 if rid in RIOS_K else .8, zorder=4)
    ax.annotate(rid, (p.x_utm.iloc[0], p.y_utm.iloc[0]), xytext=(4, 3),
                textcoords="offset points", fontsize=9,
                bbox=dict(facecolor="white", edgecolor="none", alpha=.8))
manuales_mapa = kp_litologia.loc[
    kp_litologia.ID.astype(str).isin(ids_region) & kp_litologia.tipo.eq("manual_reubicado")
].drop_duplicates(["ID", "id_kp"])
ax.scatter(manuales_mapa.x_utm, manuales_mapa.y_utm, marker="*", color="black",
           s=65, edgecolor="white", linewidth=.5, zorder=6)
ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax), xlabel="Este UTM (m)", ylabel="Norte UTM (m)",
       title="Análisis regional: cauces, selección K y KP manuales · EPSG:32718")
ax.set_aspect("equal"); ax.ticklabel_format(style="plain", useOffset=False); ax.grid(alpha=.15)
ax.legend(handles=[
    Line2D([0], [0], color="#146fa1", lw=2, label="Ríos de calibración K"),
    Line2D([0], [0], color="0.4", lw=1, label="Otros ríos regionales"),
    Line2D([0], [0], color="#bd392f", lw=1.2, label="Fallas cartografiadas"),
    Line2D([0], [0], marker="*", color="black", linestyle="none", markersize=10, label="KP manual"),
], loc="best", fontsize=9)
if GUARDAR_REGIONAL:
    sello = datetime.now(ZoneInfo("America/Los_Angeles")).strftime("%Y%m%d_%H%M%S_%f")
    carpeta = Path(SALIDA); carpeta.mkdir(parents=True, exist_ok=True)
    tabla_regional.to_csv(carpeta/f"resumen_regional_{sello}.csv", index=False, encoding="utf-8-sig")
    fig.savefig(carpeta/f"mapa_regional_{sello}.png", dpi=300, bbox_inches="tight")
    print("Guardado en:", carpeta)
plt.show(); plt.close(fig)


# Notebook cell 15
if GUARDAR_TABLAS:
    SALIDA.mkdir(parents=True, exist_ok=True)

    tablas = {
        "qc_perfiles": qc_perfiles,
        "resumen_rios": resumen_rios,
        "segmentos_litologia": pd.concat(
            segmentos.values(), ignore_index=True
        ),
        "cambios_cartograficos": cambios_geo,
        "propuestas_aluviales": propuestas_aluviales,
        "kp_litologia": kp_litologia,
        "kp_rechazados": kp_rechazados,
    }

    for nombre, tabla in tablas.items():
        tabla.to_csv(
            SALIDA / f"{nombre}_{FECHA}.csv",
            index=False,
            encoding="utf-8-sig",
        )

    ruta_propuestas = (
        SALIDA / f"propuestas_aluviales_{FECHA}.csv"
    )

    # Este JSON permite comprobar después que las propuestas
    # corresponden a los mismos perfiles.
    procedencia = {
        "fecha": FECHA,
        "theta": THETA,
        "A0_m2": A0,
        "fuentes": FUENTES,
        "perfiles_hash": {
            rid: MANIFEST_PERFILES["rios"][rid]["perfil_hash"]
            for rid in perfiles
        },
        "propuestas_path": str(ruta_propuestas),
        "propuestas_sha256": hashlib.sha256(
            ruta_propuestas.read_bytes()
        ).hexdigest(),
    }

    ruta_procedencia = SALIDA / f"procedencia_{FECHA}.json"

    ruta_procedencia.write_text(
        json.dumps(procedencia, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("Tablas guardadas en:", SALIDA)
    print("\nPara usar estas propuestas en K, copia esta ruta:")
    print(ruta_procedencia)
