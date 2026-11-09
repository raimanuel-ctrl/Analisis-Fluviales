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

ETAPA = "K_t"

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
IDS_CARGAR = ["3", "4", "7", "8", "9", "10"]

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

INCLUIR_CANDIDATOS = False
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
# K1 — FUENTES Y 30 ESCENARIOS
import re

RIOS_K = ["3", "4", "7", "8", "9", "10"]
MIN_NODOS_K = 5

if not (np.isclose(THETA, 0.40) and np.isclose(N, 1.0)):
    raise ValueError("Estos bloques requieren THETA=0.40 y N=1.")

RUTA_AGES = BASE / "Excel y csv" / "Ages_Final.xlsx"

RUTA_TERRAZAS = Path(
    r"C:\datos\Memoria_Codex"
    r"\organizacion_memoria_2_0\10_resultados\terrazas_JaraMunoz"
    r"\terrazas_170_UTM18S_con_edad_referencia_Jara_Munoz_2015.csv"
)

RUTA_OUTLETS = BASE / (
    "Productos/Redes_de_drenaje/"
    "20260929_144303_906995/rios_outlets.shp"
)

# Registrar las fuentes. El shapefile documenta la entrada confirmada
# por ti; este bloque no reconstruye ni vuelve a recortar los perfiles.
for ruta in [RUTA_AGES, RUTA_TERRAZAS, RUTA_MANUALES]:
    if not ruta.exists():
        raise FileNotFoundError(ruta)
    registrar_archivo(ruta)

if not RUTA_OUTLETS.exists():
    raise FileNotFoundError(RUTA_OUTLETS)

for extension in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
    archivo = RUTA_OUTLETS.with_suffix(extension)
    if archivo.exists():
        registrar_archivo(archivo)

# Nombres de presentación; los perfiles mantienen sus IDs originales.
NOMBRES.update({
    "3": "Estero Quebrada Honda",
    "4": "Estero Malhueco",
    "7": "Estero de Chanco",
    "8": "Río Rahue",
    "9": "Río Curanilahue",
    "10": "Río Curanipe",
})

# OD transcrito de la revisión documental previa.
# Ages_Final no contiene una columna de OD.
OD = {"CONT 11": 11.70, "PUHE 10": 14.80, "PUHE 11": 8.98}
FUENTE_OD = "Revisión documental de muestras, 20261005"

ages = pd.read_excel(RUTA_AGES, sheet_name="Sheet1")
ages.columns = ages.columns.astype(str).str.strip()

requeridas = [
    "Sample ID", "Elevation a.s.l. (m)",
    "AGE (ka)", "error", "Urate (m/ka)",
]
faltantes = set(requeridas) - set(ages.columns)
if faltantes:
    raise ValueError(f"Faltan columnas en Ages_Final: {faltantes}")

def normalizar_codigo(valor):
    return re.sub(r"[^A-Z0-9]", "", str(valor).upper())

codigos = ages["Sample ID"].map(normalizar_codigo)
MUESTRAS = {}

for nombre in OD:
    codigo = normalizar_codigo(f"JJ12-{nombre}-IRSL1")
    coincidencias = ages.loc[codigos.eq(codigo)]

    if len(coincidencias) != 1:
        raise ValueError(
            f"{nombre}: se esperaba una fila y hay {len(coincidencias)}."
        )

    r = coincidencias.iloc[0]
    datos = {
        "edad_ka": float(r["AGE (ka)"]),
        "error_edad_ka": float(r["error"]),
        "U_m_ka": float(r["Urate (m/ka)"]),
        "cota_muestra_m": float(r["Elevation a.s.l. (m)"]),
        "OD_pct": OD[nombre],
    }

    if not np.isfinite(list(datos.values())).all():
        raise ValueError(f"{nombre}: datos no finitos.")
    if (
        datos["edad_ka"] <= 0
        or datos["error_edad_ka"] < 0
        or datos["U_m_ka"] <= 0
        or not 0 <= datos["OD_pct"] < 15
    ):
        raise ValueError(f"{nombre}: revisar edad, U u OD.")

    MUESTRAS[nombre] = datos

terrazas = pd.read_csv(RUTA_TERRAZAS, encoding="utf-8-sig")
terrazas["terrace_id"] = terrazas["terrace_id"].astype(str).str.strip()

if terrazas["terrace_id"].duplicated().any():
    raise ValueError("Hay IDs de terraza duplicados en el catálogo.")

terrazas = terrazas.set_index("terrace_id")

# Selección exacta que entregaste.
REFERENCIAS = {
    "3": ["JM2015_028", "JM2015_030"],
    "4": ["JM2015_030", "JM2015_031"],
    "7": ["JM2015_045", "JM2015_157"],
    "8": ["JM2015_046", "JM2015_157"],
    "9": ["JM2015_157", "JM2015_158"],
    "10": ["JM2015_163", "JM2015_158"],
}

filas = []
for rid in RIOS_K:
    muestras_rio = (
        ["CONT 11"] if rid in ["3", "4"]
        else ["PUHE 10", "PUHE 11"]
    )

    for muestra in muestras_rio:
        datos = MUESTRAS[muestra]

        for referencia in REFERENCIAS[rid] + [muestra]:
            es_muestra = referencia == muestra

            if es_muestra:
                cota = datos["cota_muestra_m"]
                mis = ""
                fuente_cota = RUTA_AGES.name
            else:
                if referencia not in terrazas.index:
                    raise ValueError(f"No existe {referencia} en el catálogo.")
                terraza = terrazas.loc[referencia]
                cota = float(terraza["elevacion_m"])
                mis = str(terraza["mis_asignado"])
                fuente_cota = RUTA_TERRAZAS.name

            if not np.isfinite(cota):
                raise ValueError(f"Cota inválida: {referencia}")

            filas.append({
                "escenario_id": (
                    f"R{rid}_{muestra.replace(' ', '')}_"
                    f"{referencia.replace(' ', '')}"
                ),
                "ID": rid,
                "rio": NOMBRES[rid],
                "muestra": muestra,
                "referencia": referencia,
                "tipo_cota": "muestra" if es_muestra else "terraza",
                "cota_m": cota,
                **datos,
                "MIS_terraza": mis,
                "fuente_cota": fuente_cota,
            })

calibracion = pd.DataFrame(filas)

assert len(calibracion) == 30
assert calibracion.escenario_id.is_unique

# Vacío: primer cruce para cada escenario.
# Para cambiar uno, copia su escenario_id de la tabla:
# CRUCES_ELEGIDOS = {"R9_PUHE10_JM2015_158": 2}
CRUCES_ELEGIDOS = {}

display(calibracion[[
    "escenario_id", "rio", "muestra", "referencia",
    "tipo_cota", "cota_m", "edad_ka", "U_m_ka",
]])

# Notebook cell 8
# K2 — CRUCES Y AJUSTE SIN RECORTES ADICIONALES

def buscar_cruces(p, cota):
    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    puntos, planos = [], []
    tolerancia_z = 1e-8

    for i in range(len(s) - 1):
        if not np.isfinite(z[i:i+2]).all():
            continue

        d0, d1 = z[i] - cota, z[i+1] - cota
        igual0 = abs(d0) <= tolerancia_z
        igual1 = abs(d1) <= tolerancia_z

        if igual0 and igual1:
            planos.append((float(s[i]), float(s[i+1])))
            puntos.extend([s[i], s[i+1]])
        elif igual0:
            puntos.append(s[i])
        elif igual1:
            puntos.append(s[i+1])
        elif d0 * d1 < 0:
            puntos.append(
                s[i] + (cota - z[i])
                * (s[i+1] - s[i]) / (z[i+1] - z[i])
            )

    # Quitar únicamente duplicados numéricos, no cruces cercanos.
    unicos = []
    for valor in sorted(puntos):
        if not unicos or abs(valor - unicos[-1]) > 1e-6:
            unicos.append(float(valor))
    return unicos, planos


def ajustar_hasta_cruce(p, fin):
    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    chi = p.chi_040_m.to_numpy(float)

    if not np.isfinite(fin) or not 0 < fin <= s[-1]:
        raise ValueError("El cruce no define un tramo de longitud positiva.")

    hasta = min(np.searchsorted(s, fin), len(s)-1)
    if not np.isfinite(z[:hasta+1]).all():
        raise ValueError("El tramo contiene elevaciones faltantes.")

    dentro = s <= fin
    n_originales = int(dentro.sum())
    if n_originales < MIN_NODOS_K:
        raise ValueError(
            f"Solo {n_originales} nodos; mínimo {MIN_NODOS_K}."
        )

    posiciones = np.unique(np.r_[s[dentro], fin])
    x = np.interp(posiciones, s, chi)
    y = np.interp(posiciones, s, z)

    xc, yc = x - x.mean(), y - y.mean()
    denominador = float(np.sum(xc**2))
    if denominador <= 0:
        raise ValueError("No hay variación de χ.")

    pendiente = float(np.sum(xc * yc) / denominador)
    if not np.isfinite(pendiente) or pendiente <= 0:
        raise ValueError("Pendiente z–χ no positiva.")

    intercepto = float(y.mean() - pendiente * x.mean())
    residuos = y - (intercepto + pendiente * x)
    total = float(np.sum(yc**2))

    return {
        "Mchi": pendiente,
        "intercepto_m": intercepto,
        "ksn_040": pendiente * A0**THETA,
        "R2": 1 - float(np.sum(residuos**2)) / total if total > 0 else np.nan,
        "RMSE_m": float(np.sqrt(np.mean(residuos**2))),
        "n_nodos_originales": n_originales,
        "fin_m": float(fin),
    }


limites = {}
filas_cruces = []

for _, r in calibracion.iterrows():
    puntos, planos = buscar_cruces(perfiles[r.ID], r.cota_m)
    limites[r.escenario_id] = (puntos, planos)
    elegido = CRUCES_ELEGIDOS.get(r.escenario_id, 1)

    for numero, distancia in enumerate(puntos, 1):
        filas_cruces.append({
            "escenario_id": r.escenario_id,
            "ID": r.ID,
            "referencia": r.referencia,
            "cota_m": r.cota_m,
            "cruce_id": numero,
            "distancia_m": distancia,
            "n_cruces": len(puntos),
            "seleccionado": numero == elegido,
            "en_segmento_plano": any(
                a - 1e-6 <= distancia <= b + 1e-6
                for a, b in planos
            ),
        })

cruces_K = pd.DataFrame(filas_cruces, columns=[
    "escenario_id", "ID", "referencia", "cota_m",
    "cruce_id", "distancia_m", "n_cruces",
    "seleccionado", "en_segmento_plano",
])

print("Cruces registrados:", len(cruces_K))
# Descomenta para revisar la tabla completa:
# display(cruces_K)

# Notebook cell 10
# K3 — UN RESULTADO POR ESCENARIO SELECCIONADO

filas_K = []

for _, r in calibracion.iterrows():
    fila = r.to_dict()
    fila.update({
        "estado": "no_calculable",
        "detalle": "",
        "K": np.nan,
        "ksn_040": np.nan,
        "Mchi": np.nan,
        "R2": np.nan,
        "RMSE_m": np.nan,
        "fin_m": np.nan,
        "n_nodos_originales": 0,
        "unidad_K": "m^0.2/año",
    })

    puntos, planos = limites[r.escenario_id]
    numero = CRUCES_ELEGIDOS.get(r.escenario_id, 1)
    fila["cruce_id"] = numero
    fila["n_cruces"] = len(puntos)

    try:
        if not puntos:
            raise ValueError("No se encontró cruce con esta cota.")

        if not isinstance(numero, (int, np.integer)):
            raise ValueError("El número de cruce debe ser entero.")
        if not 1 <= numero <= len(puntos):
            raise ValueError("El cruce elegido no existe.")

        fin = puntos[numero - 1]

        if any(a - 1e-6 <= fin <= b + 1e-6 for a, b in planos):
            raise ValueError(
                "Cruce sobre segmento plano a la cota: revisar manualmente."
            )

        ajuste = ajustar_hasta_cruce(perfiles[r.ID], fin)
        fila.update(ajuste)
        fila["K"] = (r.U_m_ka / 1000.0) / ajuste["ksn_040"]
        fila["estado"] = "calculado_condicional"
        fila["detalle"] = (
            "Cruce único" if len(puntos) == 1
            else f"Cruce {numero} de {len(puntos)}"
        )

    except ValueError as exc:
        fila["detalle"] = str(exc)

    filas_K.append(fila)

resultados_K = pd.DataFrame(filas_K)

resumen_K = resultados_K[[
    "escenario_id", "ID", "rio", "muestra", "referencia",
    "cota_m", "cruce_id", "fin_m", "K", "estado", "detalle",
]].copy()

display(resumen_K.style.format({
    "cota_m": "{:.2f}",
    "fin_m": "{:.1f}",
    "K": "{:.3e}",
}, na_rep="—"))

print(
    "Calculables:",
    resultados_K.estado.eq("calculado_condicional").sum(),
    "de", len(resultados_K),
)

# Notebook cell 12
# K4 — SOLO OUTLET ACTUAL Y KNICKPOINTS MANUALES

if not RUTA_MANUALES.exists():
    raise FileNotFoundError(
        "Falta el archivo manual; no se puede interpretar como ausencia de KP."
    )

kp_manuales = kp.loc[
    kp.tipo.eq("manual_reubicado") & kp.ID.isin(RIOS_K)
].copy()

# Evitar repetir un mismo punto físico con dos nombres.
kp_manuales = kp_manuales.drop_duplicates(["ID", "indice"])

filas_t = []

for _, r in resultados_K.iterrows():
    base = {
        "escenario_id": r.escenario_id,
        "ID": r.ID,
        "rio": r.rio,
        "muestra": r.muestra,
        "referencia": r.referencia,
        "K": r.K,
        "origen": "outlet_recortado",
        "origen_m": 0.0,
        "id_kp": None,
        "distancia_kp_m": np.nan,
        "delta_chi_m": np.nan,
        "t_ka": np.nan,
        "fuera_tramo_calibracion": None,
    }

    if r.estado != "calculado_condicional":
        filas_t.append({**base, "estado": "K_no_calculable"})
        continue

    puntos = kp_manuales.loc[kp_manuales.ID.eq(r.ID)]
    if puntos.empty:
        hay_rechazados = (
            kp_rechazados.ID.astype(str).eq(r.ID)
            & kp_rechazados.origen.eq("manual")
        ).any()

        filas_t.append({
            **base,
            "estado": (
                "seleccion_manual_rechazada"
                if hay_rechazados else "sin_seleccion_manual"
            ),
        })
        continue

    p = perfiles[r.ID]
    chi_origen = float(p.chi_040_m.iloc[0])

    for _, punto in puntos.iterrows():
        delta_chi = float(punto.chi_040_m) - chi_origen
        if not np.isfinite(delta_chi) or delta_chi < 0:
            raise ValueError(f"χ incompatible para {punto.id_kp}")

        filas_t.append({
            **base,
            "id_kp": punto.id_kp,
            "distancia_kp_m": float(punto.distancia_m),
            "delta_chi_m": delta_chi,
            "t_ka": delta_chi / (r.K * A0**THETA) / 1000.0,
            "fuera_tramo_calibracion": bool(
                punto.distancia_m > r.fin_m
            ),
            "estado": "calculado_condicional",
        })

tiempos = pd.DataFrame(filas_t)

display(tiempos[[
    "escenario_id", "id_kp", "K", "t_ka",
    "fuera_tramo_calibracion", "estado",
]].style.format({
    "K": "{:.3e}",
    "t_ka": "{:.2f}",
}, na_rep="—"))

# Notebook cell 14
# K5 — EXPORTACIÓN Y FUNCIONES DE VISUALIZACIÓN
import matplotlib.pyplot as plt

if GUARDAR_TABLAS:
    SALIDA.mkdir(parents=True, exist_ok=True)

    tablas = {
        "calibracion_30_escenarios": calibracion,
        "resultados_K": resultados_K,
        "resumen_K": resumen_K,
        "cruces_K": cruces_K,
        "tiempos": tiempos,
        "kp_manuales": kp_manuales,
        "kp_rechazados": kp_rechazados,
    }

    for nombre, tabla in tablas.items():
        tabla.to_csv(
            SALIDA / f"{nombre}_{FECHA}.csv",
            index=False,
            encoding="utf-8-sig",
            float_format="%.12g",
        )

    parametros = {
        "fecha": FECHA,
        "theta": THETA,
        "A0_m2": A0,
        "n": N,
        "min_nodos_K": MIN_NODOS_K,
        "metodo_K": "(U_m_ka / 1000) / ksn",
        "unidad_K": "m^0.2/año",
        "origen_t": "outlet_recortado",
        "recortes_adicionales": False,
        "perfiles_ya_recortados": "Confirmado por el usuario",
        "outlets_declarados": str(RUTA_OUTLETS),
        "cruces_elegidos": CRUCES_ELEGIDOS,
        "fuente_OD": FUENTE_OD,
        "fuentes_sha256": FUENTES,
        "escenarios": calibracion.to_dict(orient="records"),
    }

    (SALIDA / f"parametros_{FECHA}.json").write_text(
        json.dumps(
            parametros, indent=2, ensure_ascii=False, allow_nan=False
        ),
        encoding="utf-8",
    )
    print("Resultados guardados en:", SALIDA)


def seleccionar_rios(rios=None):
    """Vacío muestra los seis; admite un ID o una colección de IDs."""
    if rios is None:
        return RIOS_K.copy()
    if isinstance(rios, (str, int)):
        rios = [rios]
    elegidos = {str(r) for r in rios}
    if not elegidos:
        return RIOS_K.copy()
    if elegidos - set(RIOS_K):
        raise ValueError("Selecciona solo IDs 3, 4, 7, 8, 9 y 10.")
    return [rid for rid in RIOS_K if rid in elegidos]


def mostrar_K(rios=None):
    tabla = resumen_K.loc[
        resumen_K.ID.isin(seleccionar_rios(rios))
    ]
    display(tabla.style.format({
        "cota_m": "{:.2f}",
        "fin_m": "{:.1f}",
        "K": "{:.3e}",
    }, na_rep="—"))


def mostrar_t(rios=None):
    tabla = tiempos.loc[
        tiempos.ID.isin(seleccionar_rios(rios)),
        ["escenario_id", "id_kp", "K", "t_ka",
         "fuera_tramo_calibracion", "estado"],
    ]
    display(tabla.style.format({
        "K": "{:.3e}", "t_ka": "{:.2f}",
    }, na_rep="—"))


def mostrar_cruces(rios=None, perfil_completo=False, guardar=False):
    ids = seleccionar_rios(rios)
    columnas = min(3, len(ids))
    filas = int(np.ceil(len(ids) / columnas))

    fig, axes = plt.subplots(
        filas, columnas,
        figsize=(6 * columnas, 4.5 * filas),
        squeeze=False,
        constrained_layout=True,
    )

    for ax, rid in zip(axes.flat, ids):
        p = perfiles[rid]
        ax.plot(
            p.distancia_m / 1000, p.elevacion_m,
            color="#24485C", lw=1.5, label="Cauce",
        )

        referencias = (
            calibracion.loc[calibracion.ID.eq(rid)]
            .drop_duplicates("referencia")
        )
        colores = plt.get_cmap("tab10")
        distancias_seleccionadas = []

        for j, (_, ref) in enumerate(referencias.iterrows()):
            color = colores(j % 10)
            etiqueta = (
                f"{ref.tipo_cota}: {ref.referencia} "
                f"({ref.cota_m:.2f} m)"
            )
            ax.axhline(
                ref.cota_m, color=color, lw=1.2,
                ls="--", label=etiqueta,
            )

            puntos, _ = limites[ref.escenario_id]
            ax.scatter(
                np.asarray(puntos) / 1000,
                np.full(len(puntos), ref.cota_m),
                facecolors="none", edgecolors=[color],
                s=32, zorder=4,
            )

            # Ambas muestras pueden compartir una referencia.
            # Si sus cruces elegidos difieren, se muestran ambos.
            elegidos = resultados_K.loc[
                resultados_K.ID.eq(rid)
                & resultados_K.referencia.eq(ref.referencia)
                & resultados_K.estado.eq("calculado_condicional"),
                "fin_m",
            ].dropna().unique()

            distancias_seleccionadas.extend(elegidos)
            ax.scatter(
                elegidos / 1000,
                np.full(len(elegidos), ref.cota_m),
                marker="x", color="black", s=55, zorder=5,
            )

        manuales = kp_manuales.loc[kp_manuales.ID.eq(rid)]
        if not manuales.empty:
            ax.scatter(
                manuales.distancia_m / 1000,
                manuales.elevacion_m,
                marker="*", color="black", s=85,
                label="KP manual", zorder=6,
            )

        if not perfil_completo and distancias_seleccionadas:
            fin_vista = min(
                float(p.distancia_m.iloc[-1]),
                1.2 * max(distancias_seleccionadas),
            )
            ax.set_xlim(0, fin_vista / 1000)

            local = p.loc[p.distancia_m <= fin_vista, "elevacion_m"]
            alturas = np.r_[
                local.to_numpy(float), referencias.cota_m.to_numpy(float)
            ]
            margen = max(5.0, 0.1 * np.ptp(alturas))
            ax.set_ylim(np.min(alturas)-margen, np.max(alturas)+margen)

        ax.set(
            title=f"{rid} · {NOMBRES[rid]}",
            xlabel="Distancia desde outlet recortado (km)",
            ylabel="Elevación (m)",
        )
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7)

    for ax in axes.flat[len(ids):]:
        ax.set_visible(False)

    if guardar:
        SALIDA.mkdir(parents=True, exist_ok=True)
        etiqueta = "_".join(ids)
        vista = "completo" if perfil_completo else "calibracion"
        fig.savefig(
            SALIDA / f"cruces_{etiqueta}_{vista}_{FECHA}.png",
            dpi=300,
        )

    plt.show()
    plt.close(fig)

print("Funciones disponibles: mostrar_cruces(), mostrar_K(), mostrar_t().")

# Notebook cell 15
# Los seis ríos, aproximadamente seis gráficos por página:
mostrar_cruces()

# Solo los ríos indicados:
# mostrar_cruces({"3", "4"})
# mostrar_cruces({"9"})

# Perfil completo para ver también los KP más alejados:
# mostrar_cruces({"9"}, perfil_completo=True)

# Guardar la figura con fecha:
# mostrar_cruces(guardar=True)

# Tablas:

#mostrar_K()
# mostrar_K({"9"})

#mostrar_t({"9"})

# Notebook cell 17
# K6 — DIAGNÓSTICO POR RÍO
# Requiere haber ejecutado K1–K5.
# No modifica resultados_K, perfiles ni la selección de cruces.

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def diagnostico_calibracion(rid, guardar=False):
    rid = str(rid)

    if rid not in RIOS_K:
        raise ValueError("Selecciona un río: 3, 4, 7, 8, 9 o 10.")

    p = perfiles[rid]
    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    chi = p.chi_040_m.to_numpy(float)

    resultados_rio = resultados_K.loc[
        resultados_K.ID.eq(rid)
    ].copy()

    validos = resultados_rio.loc[
        resultados_rio.estado.eq("calculado_condicional")
    ].copy()

    no_calculables = resultados_rio.loc[
        ~resultados_rio.estado.eq("calculado_condicional"),
        ["escenario_id", "estado", "detalle"],
    ]

    if not no_calculables.empty:
        print("Escenarios que no entran al diagnóstico:")
        display(no_calculables)

    if validos.empty:
        print("Este río no tiene calibraciones calculables.")
        return

    # Un grupo por geometría.
    # PUHE 10 y PUHE 11 comparten el ajuste cuando la cota
    # y la distancia del cruce son exactamente las mismas.
    grupos = list(
        validos.groupby(["cota_m", "fin_m"], sort=True)
    )

    manuales = kp_manuales.loc[kp_manuales.ID.eq(rid)].copy()

    fig, axes = plt.subplots(
        1, 3, figsize=(18, 5), constrained_layout=True
    )
    ax_completo, ax_ajuste, ax_residuos = axes

    ax_completo.plot(
        chi, z, color="0.65", lw=1.2,
        label="Perfil completo", zorder=1,
    )

    colores = plt.get_cmap("tab10")
    filas_diagnostico = []
    leyenda_tramos = []
    limites_chi = []

    for numero, ((cota, fin), grupo) in enumerate(grupos, 1):
        color = colores((numero - 1) % 10)

        # Reconstruir los mismos puntos usados por el ajuste K2:
        # nodos originales + extremo interpolado.
        dentro = s <= fin
        posiciones = np.unique(np.r_[s[dentro], fin])
        x = np.interp(posiciones, s, chi)
        y = np.interp(posiciones, s, z)

        ajuste = ajustar_hasta_cruce(p, float(fin))
        prediccion = (
            ajuste["intercepto_m"] + ajuste["Mchi"] * x
        )
        residuos = y - prediccion

        referencias = " / ".join(
            sorted(grupo.referencia.astype(str).unique())
        )
        muestras = " / ".join(
            sorted(grupo.muestra.astype(str).unique())
        )
        etiqueta = f"T{numero}: {referencias} · {cota:.2f} m"
        limites_chi.append(float(x[-1]))

        # Panel 1: recta únicamente dentro del tramo calibrado.
        # No dibujamos una extrapolación como si fuera un ajuste.
        ax_completo.plot(
            x, prediccion, color=color, lw=2.2, zorder=3
        )
        ax_completo.scatter(
            [x[-1]], [cota], color=color, s=38, zorder=4
        )

        # Panel 2: datos y recta en el tramo inferior.
        ax_ajuste.plot(
            x, y, color=color, lw=0.9, alpha=0.65
        )
        ax_ajuste.plot(
            x, prediccion, color=color, lw=2, ls="--"
        )
        ax_ajuste.scatter(
            [x[-1]], [cota], color=color, s=35, zorder=4
        )

        # Panel 3: residuos sin suavizar.
        ax_residuos.plot(
            x, residuos, color=color, lw=1.1
        )

        leyenda_tramos.append(
            Line2D([0], [0], color=color, lw=2, label=etiqueta)
        )

        filas_diagnostico.append({
            "tramo": f"T{numero}",
            "referencias": referencias,
            "muestras": muestras,
            "cota_m": float(cota),
            "longitud_m": float(fin),
            "n_nodos": ajuste["n_nodos_originales"],
            "ksn_040": ajuste["ksn_040"],
            "R2": ajuste["R2"],
            "RMSE_m": ajuste["RMSE_m"],
            "residuo_max_abs_m": float(
                np.max(np.abs(residuos))
            ),
            # Extensión relativa del tramo respecto del perfil.
            "fraccion_longitud_perfil": float(fin / s[-1]),
        })

    # KP manuales sobre el perfil completo.
    if not manuales.empty:
        ax_completo.scatter(
            manuales.chi_040_m,
            manuales.elevacion_m,
            marker="*", s=130, color="black",
            label="KP manual", zorder=6,
        )

        for _, punto in manuales.iterrows():
            ax_completo.annotate(
                str(punto.id_kp),
                (punto.chi_040_m, punto.elevacion_m),
                xytext=(5, 7),
                textcoords="offset points",
                fontsize=7,
            )
    else:
        ax_completo.text(
            0.03, 0.97, "Sin KP manual identificado",
            transform=ax_completo.transAxes,
            va="top", fontsize=9,
        )

    # Delimitar visualmente el sector ampliado en el panel 2.
    ax_completo.axvspan(
        0, max(limites_chi), color="0.8", alpha=0.15,
        label="Sector mostrado en detalle",
    )

    ax_completo.set(
        title="1. Contexto del perfil completo",
        xlabel="χ (m)",
        ylabel="Elevación (m)",
    )
    ax_completo.legend(fontsize=7, loc="best")

    ax_ajuste.set(
        title="2. Tramos y ajustes",
        xlabel="χ (m)",
        ylabel="Elevación (m)",
    )
    ax_ajuste.legend(
        handles=leyenda_tramos, fontsize=7, loc="best"
    )

    ax_residuos.axhline(
        0, color="black", lw=1, ls="--"
    )
    ax_residuos.set(
        title="3. Residuos del ajuste",
        xlabel="χ (m)",
        ylabel="Elevación observada − ajustada (m)",
    )

    for ax in axes:
        ax.grid(alpha=0.2)

    fig.suptitle(
        f"{rid} · {NOMBRES[rid]} — diagnóstico de calibración\n"
        "Panel 2: líneas continuas = datos; discontinuas = ajustes",
        fontsize=12,
    )

    diagnostico = pd.DataFrame(filas_diagnostico)

    # Tabla que hace explícita la extrapolación espacial hasta cada KP.
    extrapolaciones = []
    for _, tramo in diagnostico.iterrows():
        for _, punto in manuales.iterrows():
            distancia_kp = float(punto.distancia_m)
            extrapolaciones.append({
                "tramo": tramo.tramo,
                "id_kp": punto.id_kp,
                "fin_calibracion_m": tramo.longitud_m,
                "distancia_kp_m": distancia_kp,
                "recorrido_fuera_calibracion_m": max(
                    0.0, distancia_kp - tramo.longitud_m
                ),
            })

    extrapolacion = pd.DataFrame(extrapolaciones)

    if guardar:
        SALIDA.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            SALIDA / f"diagnostico_rio_{rid}_{FECHA}.png",
            dpi=300,
        )
        diagnostico.to_csv(
            SALIDA / f"diagnostico_rio_{rid}_{FECHA}.csv",
            index=False, encoding="utf-8-sig",
        )
        if not extrapolacion.empty:
            extrapolacion.to_csv(
                SALIDA / f"extrapolacion_rio_{rid}_{FECHA}.csv",
                index=False, encoding="utf-8-sig",
            )

    plt.show()
    plt.close(fig)

    print("Diagnóstico por tramo geométrico:")
    display(diagnostico.style.format({
        "cota_m": "{:.2f}",
        "longitud_m": "{:.1f}",
        "ksn_040": "{:.2f}",
        "R2": "{:.3f}",
        "RMSE_m": "{:.2f}",
        "residuo_max_abs_m": "{:.2f}",
        "fraccion_longitud_perfil": "{:.1%}",
    }))

    if not extrapolacion.empty:
        print("Recorrido hasta el KP fuera del tramo calibrado:")
        display(extrapolacion.style.format({
            "fin_calibracion_m": "{:.1f}",
            "distancia_kp_m": "{:.1f}",
            "recorrido_fuera_calibracion_m": "{:.1f}",
        }))


# COMENZAMOS CON HONDA:
diagnostico_calibracion("3")
# Luego, ejecutar uno por uno:
# diagnostico_calibracion("8")
# diagnostico_calibracion("9")

# Guardar figura y tablas con fecha:
# diagnostico_calibracion("3", guardar=True)

# Notebook cell 19
# K7 — LITOLOGÍA, CALIBRACIÓN Y RECORRIDO HASTA EL KP
# Requiere K1–K6. No cambia K ni vuelve a recortar los perfiles.

from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.colors import is_color_like

FECHA_LITO = "20261006_132435_105931"
CARPETA_LITO = BASE / "Productos" / "litologia" / FECHA_LITO

RUTA_SEGMENTOS = (
    CARPETA_LITO / f"segmentos_litologia_{FECHA_LITO}.csv"
)
RUTA_PROCEDENCIA = (
    CARPETA_LITO / f"procedencia_{FECHA_LITO}.json"
)

for ruta in [RUTA_SEGMENTOS, RUTA_PROCEDENCIA]:
    if not ruta.exists():
        raise FileNotFoundError(ruta)
    registrar_archivo(ruta)

procedencia = json.loads(
    RUTA_PROCEDENCIA.read_text(encoding="utf-8-sig")
)

segmentos_geo = pd.read_csv(
    RUTA_SEGMENTOS,
    dtype={"ID": str},
    encoding="utf-8-sig",
    float_precision="round_trip",
)

# Comprobar la procedencia y la geometría:
# no mezclar litología exportada sobre otro perfil.
for rid in RIOS_K:
    esperado = procedencia["perfiles_hash"].get(rid)
    actual = MANIFEST_PERFILES["rios"][rid]["perfil_hash"]

    if esperado != actual:
        raise ValueError(f"Río {rid}: litología de otro perfil.")

    q = segmentos_geo.loc[
        segmentos_geo.ID.eq(rid)
    ].sort_values("s0")
    p = perfiles[rid]

    if len(q) != len(p) - 1:
        raise ValueError(f"Río {rid}: número de segmentos incompatible.")

    for campo_geo, campo_perfil, extremo in [
        ("s0", "distancia_m", "inicio"),
        ("s1", "distancia_m", "fin"),
        ("z0", "elevacion_m", "inicio"),
        ("z1", "elevacion_m", "fin"),
    ]:
        valores = p[campo_perfil].to_numpy(float)
        valores = valores[:-1] if extremo == "inicio" else valores[1:]

        if not np.allclose(
            q[campo_geo].to_numpy(float), valores,
            rtol=0, atol=1e-6, equal_nan=True,
        ):
            raise ValueError(
                f"Río {rid}: diferencias en {campo_geo}."
            )


def revisar_litologia(rid, guardar=False):
    rid = str(rid)
    if rid not in RIOS_K:
        raise ValueError("Usa un ID entre 3, 4, 7, 8, 9 y 10.")

    p = perfiles[rid]
    q = (
        segmentos_geo.loc[segmentos_geo.ID.eq(rid)]
        .sort_values("s0").copy()
    )
    manuales = kp_manuales.loc[kp_manuales.ID.eq(rid)].copy()

    validos = resultados_K.loc[
        resultados_K.ID.eq(rid)
        & resultados_K.estado.eq("calculado_condicional")
    ].copy()

    if validos.empty:
        print("No hay calibraciones válidas para este río.")
        return

    # Agrupar escenarios con el mismo límite.
    # Una terraza usada con ambas PUHE no se dibuja dos veces.
    limites_vista = []
    for numero, (fin, grupo) in enumerate(
        validos.groupby("fin_m", sort=True), 1
    ):
        limites_vista.append({
            "tramo": f"L{numero}",
            "fin_m": float(fin),
            "referencias": " / ".join(
                sorted(grupo.referencia.unique())
            ),
        })

    # Mostrar hasta el KP más lejano o hasta el último límite,
    # con un pequeño margen. Si no hay KP, mostrar hasta ese límite.
    extremo = max(
        validos.fin_m.max(),
        manuales.distancia_m.max() if not manuales.empty else 0,
    )
    fin_vista = min(
        float(p.distancia_m.iloc[-1]), extremo * 1.05
    )

    # Recortar solo la VISUALIZACIÓN al extremo elegido.
    vista = q.loc[q.s0 < fin_vista].copy()
    vista["s1"] = np.minimum(vista.s1, fin_vista)

    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    chi = p.chi_040_m.to_numpy(float)

    vista["z0"] = np.interp(vista.s0, s, z)
    vista["z1"] = np.interp(vista.s1, s, z)
    vista["xchi0"] = np.interp(vista.s0, s, chi)
    vista["xchi1"] = np.interp(vista.s1, s, chi)

    def etiqueta_unidad(r):
        if r.estado_geo != "asignado":
            return str(r.estado_geo)
        return f"{r.codigo} · {r.unidad} · {r.fuente}"

    vista["descripcion"] = vista.apply(etiqueta_unidad, axis=1)
    colores = [
        c if isinstance(c, str) and is_color_like(c)
        else "#B7B7B7"
        for c in vista.color_hex
    ]

    fig, axes = plt.subplots(
        2, 1, figsize=(13, 9), constrained_layout=True
    )

    for ax, x0, x1, xlabel in [
        (
            axes[0],
            vista.s0.to_numpy() / 1000,
            vista.s1.to_numpy() / 1000,
            "Distancia desde outlet recortado (km)",
        ),
        (
            axes[1],
            vista.xchi0.to_numpy(),
            vista.xchi1.to_numpy(),
            "χ (m)",
        ),
    ]:
        segmentos_linea = np.stack([
            np.column_stack([x0, vista.z0]),
            np.column_stack([x1, vista.z1]),
        ], axis=1)

        # Borde oscuro para distinguir colores claros.
        ax.add_collection(LineCollection(
            segmentos_linea, colors="#333333", linewidths=4.5
        ))
        ax.add_collection(LineCollection(
            segmentos_linea, colors=colores, linewidths=3
        ))
        ax.autoscale_view()
        ax.set(xlabel=xlabel, ylabel="Elevación (m)")
        ax.grid(alpha=0.2)

        for limite in limites_vista:
            posicion = (
                limite["fin_m"] / 1000
                if ax is axes[0]
                else np.interp(limite["fin_m"], s, chi)
            )
            ax.axvline(
                posicion, color="0.35", ls="--", lw=0.9
            )
            ax.text(
                posicion, 0.98, limite["tramo"],
                transform=ax.get_xaxis_transform(),
                rotation=90, va="top", ha="right", fontsize=8,
            )

        if not manuales.empty:
            xkp = (
                manuales.distancia_m / 1000
                if ax is axes[0] else manuales.chi_040_m
            )
            ax.scatter(
                xkp, manuales.elevacion_m,
                marker="*", s=120, color="black", zorder=6,
            )

    leyenda = []
    for _, r in vista.drop_duplicates("descripcion").iterrows():
        color = (
            r.color_hex if is_color_like(r.color_hex)
            else "#B7B7B7"
        )
        leyenda.append(Line2D(
            [0], [0], color=color, lw=4, label=r.descripcion
        ))

    if not manuales.empty:
        leyenda.append(Line2D(
            [0], [0], marker="*", color="black",
            linestyle="none", markersize=10, label="KP manual",
        ))

    axes[1].legend(
        handles=leyenda, fontsize=8,
        loc="upper left", bbox_to_anchor=(0, -0.18),
    )
    fig.suptitle(f"{rid} · {NOMBRES[rid]} — litología y calibración")

    # Longitud de cada unidad por intervalo.
    # Incluye explícitamente zonas sin cobertura y solapes.
    resumen = []

    def resumir_intervalo(tramo, sector, inicio, fin, id_kp=""):
        if fin <= inicio:
            return

        solape = np.maximum(
            0,
            np.minimum(q.s1.to_numpy(), fin)
            - np.maximum(q.s0.to_numpy(), inicio),
        )
        parte = q.loc[solape > 0].copy()
        parte["longitud_m"] = solape[solape > 0]
        parte["descripcion"] = parte.apply(etiqueta_unidad, axis=1)

        for descripcion, grupo in parte.groupby("descripcion"):
            longitud = float(grupo.longitud_m.sum())
            resumen.append({
                "tramo": tramo,
                "sector": sector,
                "id_kp": id_kp,
                "unidad_y_fuente": descripcion,
                "longitud_m": longitud,
                "porcentaje_intervalo": 100 * longitud / (fin - inicio),
            })

    for limite in limites_vista:
        fin = limite["fin_m"]
        resumir_intervalo(
            limite["tramo"], "Calibración", 0.0, fin
        )
        for _, punto in manuales.iterrows():
            if punto.distancia_m > fin:
                resumir_intervalo(
                    limite["tramo"], "Fuera de calibración",
                    fin, float(punto.distancia_m), punto.id_kp,
                )

    tabla_litologia = pd.DataFrame(resumen)

    if guardar:
        SALIDA.mkdir(parents=True, exist_ok=True)
        fig.savefig(
            SALIDA / f"litologia_calibracion_rio_{rid}_{FECHA}.png",
            dpi=300, bbox_inches="tight",
        )
        tabla_litologia.to_csv(
            SALIDA / f"litologia_calibracion_rio_{rid}_{FECHA}.csv",
            index=False, encoding="utf-8-sig",
        )

    plt.show()
    plt.close(fig)

    print("Correspondencia de límites:")
    display(pd.DataFrame(limites_vista))

    print("Unidades cartográficas por sector:")
    display(tabla_litologia.style.format({
        "longitud_m": "{:.1f}",
        "porcentaje_intervalo": "{:.1f}",
    }))

    if not manuales.empty:
        print(
            "Los porcentajes describen longitud cartografiada; "
            "no prueban que K sea constante."
        )
    else:
        print("Sin KP manual: se revisa únicamente la calibración.")


# Empieza por Curanilahue:
revisar_litologia("9")

# Después, uno por uno:
revisar_litologia("8")
revisar_litologia("10")
revisar_litologia("4")
revisar_litologia("3")
revisar_litologia("7")

# Para guardar figura y tabla:
# revisar_litologia("9", guardar=True)

# Notebook cell 20
revisar_litologia("9", guardar=True)

# Notebook cell 21
# K8 — MAPA, RIBERAS Y SENSIBILIDAD A LOS CRUCES
# Reemplazar la celda K8 completa. Requiere las celdas anteriores K1–K7.
# K8 — MAPA DE CALIBRACIÓN Y ENTORNO DEL KP
import geopandas as gpd
from shapely.geometry import Point, box
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import is_color_like
import matplotlib.pyplot as plt

GEO_DIR = Path(
    r"C:\datos\Memoria_Codex\outputs"
    r"\Geologia_integrada_Chanco_Jorquera_Gajardo"
)

RUTA_GEO = GEO_DIR / "Geologia_unidades.shp"

for extension in [".shp", ".shx", ".dbf", ".prj", ".cpg"]:
    archivo = RUTA_GEO.with_suffix(extension)
    if archivo.exists():
        registrar_archivo(archivo)

geologia_mapa = gpd.read_file(RUTA_GEO)
if geologia_mapa.crs is None:
    raise ValueError("La geología no tiene CRS definido.")

geologia_mapa = geologia_mapa.to_crs("EPSG:32718")


# Capa compartida por el mapa y la comprobación bilateral.
geo = geologia_mapa
CAMPOS_GEO = ["codigo", "unidad", "grupo_lito", "fuente", "color_hex", "uso_metr", "escala"]
if set(CAMPOS_GEO) - set(geo.columns):
    raise ValueError("Faltan campos geológicos necesarios para la revisión bilateral.")

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


def consultar_litologia(x, y, rid, posicion_m):
    puntos = gpd.GeoDataFrame(geometry=[Point(float(x), float(y))], crs="EPSG:32718")
    tabla = asignar_unidades(puntos, ids=str(rid), posiciones_m=[posicion_m]).reset_index(drop=True)
    return tabla[["codigo", "unidad", "fuente", "estado_geo", "distancia_poligono_m",
                  "metodo_asignacion", "distancia_izquierda_m", "distancia_derecha_m",
                  "motivo_asignacion"]]

def _puntos_referencia_mapa_k8(seleccion, ver_terrazas=True, ver_dataciones=True):
    """Ubicaciones reales de los catálogos; no proyectar estos puntos al río."""
    import re
    registros = []

    if ver_terrazas:
        catalogo = terrazas.copy()
        if "terrace_id" not in catalogo.columns:
            catalogo = catalogo.reset_index()
        requeridas = {"terrace_id", "este_utm18s_m", "norte_utm18s_m", "elevacion_m"}
        if requeridas - set(catalogo.columns):
            raise ValueError("El catálogo de terrazas no tiene sus coordenadas UTM 18S.")
        catalogo["terrace_id"] = catalogo.terrace_id.astype(str).str.strip()
        referencias = seleccion.loc[seleccion.tipo_cota.eq("terraza"), "referencia"].unique()
        for referencia in referencias:
            filas = catalogo.loc[catalogo.terrace_id.eq(str(referencia))]
            if len(filas) != 1:
                raise ValueError(f"Terraza {referencia}: se esperaba una ubicación única.")
            r = filas.iloc[0]
            if "crs_utm" in r.index and not str(r.crs_utm).startswith("EPSG:32718"):
                raise ValueError(f"Verifica el CRS de la terraza {referencia}.")
            registros.append(dict(
                nombre=str(referencia), tipo="Terraza",
                este_m=float(r.este_utm18s_m), norte_m=float(r.norte_utm18s_m),
                cota_m=float(r.elevacion_m), fuente="Catálogo de terrazas Jara-Muñoz",
            ))

    if ver_dataciones:
        tabla = ages.copy()
        tabla.columns = tabla.columns.astype(str).str.strip()
        # Encabezados invertidos en el Ages_Final verificado:
        # Ny contiene Este (~721000); Ex contiene Norte (~6035000).
        col_este = globals().get("COLUMNA_ESTE_DATACION_K8", "Ny")
        col_norte = globals().get("COLUMNA_NORTE_DATACION_K8", "Ex")
        requeridas = {"Sample ID", "Elevation a.s.l. (m)", col_este, col_norte}
        if requeridas - set(tabla.columns):
            raise ValueError("Ages_Final no tiene las columnas de coordenadas indicadas.")
        normalizar = lambda v: re.sub(r"[^A-Z0-9]", "", str(v).upper())
        codigos = tabla["Sample ID"].map(normalizar)
        for nombre in seleccion.muestra.unique():
            codigo = normalizar(f"JJ12-{nombre}-IRSL1")
            filas = tabla.loc[codigos.eq(codigo)]
            if len(filas) != 1:
                raise ValueError(f"Muestra {nombre}: se esperaba una ubicación única.")
            r = filas.iloc[0]
            registros.append(dict(
                nombre=str(nombre), tipo="Datación",
                este_m=float(r[col_este]), norte_m=float(r[col_norte]),
                cota_m=float(r["Elevation a.s.l. (m)"]), fuente="Ages_Final.xlsx",
            ))

    puntos = pd.DataFrame(registros, columns=[
        "nombre", "tipo", "este_m", "norte_m", "cota_m", "fuente",
    ])
    if not puntos.empty:
        xy = puntos[["este_m", "norte_m"]].to_numpy(float)
        if not np.isfinite(xy).all():
            raise ValueError("Hay coordenadas de referencia no finitas.")
        if not (puntos.este_m.between(100000, 900000).all()
                and puntos.norte_m.between(0, 10000000).all()):
            raise ValueError("Coordenadas fuera de rango UTM: revisa las columnas Este/Norte.")
    return puntos


def mapa_revision(
    rid="9", radio_kp_m=1000, margen_calibracion_m=500, guardar=False,
    mostrar_rio_completo=True, mostrar_detalles=True,
    mostrar_terrazas=True, mostrar_dataciones=True, mostrar_litologia=True,
    referencias=None, muestras=None, limites_modelos=None,
):
    """Mapa general + detalles; permite límites propios de los modelos visibles."""
    from datetime import datetime
    from zoneinfo import ZoneInfo
    import textwrap

    rid = str(rid)
    if rid not in perfiles or rid not in list(map(str, RIOS_K)):
        raise ValueError("Selecciona uno de los seis ríos cargados.")
    if radio_kp_m <= 0 or margen_calibracion_m <= 0:
        raise ValueError("Los márgenes deben ser positivos.")
    if not mostrar_rio_completo and not mostrar_detalles:
        return

    p = perfiles[rid]
    s = p.distancia_m.to_numpy(float)
    este = p.x_utm.to_numpy(float)
    norte = p.y_utm.to_numpy(float)
    manuales = kp_manuales.loc[kp_manuales.ID.astype(str).eq(rid)].copy()
    seleccion = calibracion.loc[calibracion.ID.astype(str).eq(rid)].copy()
    for columna, filtro in [("referencia", referencias), ("muestra", muestras)]:
        if filtro is not None:
            filtro = [filtro] if isinstance(filtro, str) else list(filtro)
            if filtro:
                seleccion = seleccion.loc[seleccion[columna].isin(filtro)]
    if seleccion.empty:
        raise ValueError("No hay referencias para esta selección del mapa.")

    puntos = _puntos_referencia_mapa_k8(
        seleccion, mostrar_terrazas, mostrar_dataciones,
    )
    if limites_modelos is None:
        validos = resultados_K.loc[
            resultados_K.ID.astype(str).eq(rid)
            & resultados_K.estado.eq("calculado_condicional")
            & resultados_K.escenario_id.isin(seleccion.escenario_id)
        ].copy()
    else:
        # Primer cruce propio de cada modelo visible, incluso si su
        # ajuste queda pendiente. La leyenda identifica raw/suavizado.
        validos = limites_modelos.copy()

    # L1, L2... son límites sobre el río, ordenados desde el outlet.
    limites = []
    for numero, (fin, grupo) in enumerate(validos.groupby("fin_m", sort=True), 1):
        fin = float(fin)
        if not s[0] <= fin <= s[-1]:
            raise ValueError("Un límite de calibración está fuera del perfil cargado.")
        nombres = " / ".join(sorted(grupo.referencia.astype(str).unique()))
        cotas = " / ".join(f"{v:.2f}" for v in sorted(grupo.cota_m.unique()))
        limites.append(dict(
            limite=f"L{numero}", distancia_m=fin, referencias=nombres,
            cotas_m=cotas, este_m=float(np.interp(fin, s, este)),
            norte_m=float(np.interp(fin, s, norte)),
        ))
    tabla_limites = pd.DataFrame(limites, columns=[
        "limite", "distancia_m", "referencias", "cotas_m", "este_m", "norte_m",
    ])

    fin_max = None if tabla_limites.empty else float(tabla_limites.distancia_m.max())
    if fin_max is not None:
        pos_cal = np.unique(np.r_[s[s <= fin_max], fin_max])
        xcal, ycal = np.interp(pos_cal, s, este), np.interp(pos_cal, s, norte)
    else:
        xcal, ycal = este[:1], norte[:1]

    sello = datetime.now(ZoneInfo("America/Los_Angeles")).strftime("%Y%m%d_%H%M%S_%f")
    guardadas = []

    def dibujar_lamina(cajas, etiqueta_archivo):
        ancho = 8.5 if len(cajas) == 1 else 7 * len(cajas)
        fig, axes = plt.subplots(1, len(cajas), figsize=(ancho, 6.5), squeeze=False)
        unidades = {}
        referencias_visibles = set()

        for ax, (titulo, bounds) in zip(axes.flat, cajas):
            xmin, ymin, xmax, ymax = bounds
            ventana = box(xmin, ymin, xmax, ymax)
            if mostrar_litologia and rid not in RIOS_SIN_CARTA_CONFIRMADA:
                g = geologia_mapa.loc[geologia_mapa.geometry.intersects(ventana)].copy()
                if not g.empty:
                    g.geometry = g.geometry.intersection(ventana)
                    g = g.loc[~g.geometry.is_empty].copy()
                    colores = []
                    for _, u in g.iterrows():
                        color = u.color_hex
                        if not isinstance(color, str) or not is_color_like(color):
                            color = "#B7B7B7"
                        colores.append(color)
                        # Código y fuente en el mapa; nombre completo en los perfiles.
                        unidades[f"{u.codigo} · {u.fuente}"] = color
                    if not g.empty:
                        g.plot(ax=ax, color=colores, edgecolor="#777777",
                               linewidth=.5, alpha=.65, zorder=1)

            ax.plot(este, norte, color="white", lw=4, zorder=3)
            ax.plot(este, norte, color="#124B75", lw=1.8, zorder=4)
            if fin_max is not None:
                ax.plot(xcal, ycal, color="#00A6B2", lw=2.8, zorder=5)
            ax.scatter(este[0], norte[0], marker="s", color="black", s=50, zorder=8)

            for numero_punto, (_, punto) in enumerate(puntos.iterrows()):
                x, y = punto.este_m, punto.norte_m
                if not (xmin <= x <= xmax and ymin <= y <= ymax):
                    continue
                referencias_visibles.add((punto.tipo, punto.nombre))
                es_terraza = punto.tipo == "Terraza"
                color = "#7943A2" if es_terraza else "#E57A16"
                ax.scatter(x, y, marker="D" if es_terraza else "^", color=color,
                           edgecolor="white", linewidth=.7, s=80, zorder=9)
                desplazamiento = (8, 12) if numero_punto % 2 == 0 else (-8, -16)
                ax.annotate(punto.nombre, (x, y), xytext=desplazamiento,
                            ha="left" if numero_punto % 2 == 0 else "right",
                            textcoords="offset points", fontsize=8, color=color,
                            arrowprops=dict(arrowstyle="-", color=color, lw=.6),
                            bbox=dict(facecolor="white", alpha=.85, edgecolor="none"))

            for _, limite in tabla_limites.iterrows():
                x, y = limite.este_m, limite.norte_m
                if xmin <= x <= xmax and ymin <= y <= ymax:
                    ax.scatter(x, y, marker="x", color="#B52A35", s=55, lw=1.8, zorder=10)
                    ax.annotate(limite.limite, (x, y), xytext=(6, 6),
                                textcoords="offset points", fontsize=8,
                                bbox=dict(facecolor="white", alpha=.85, edgecolor="none"))
            for _, punto in manuales.iterrows():
                if xmin <= punto.x_utm <= xmax and ymin <= punto.y_utm <= ymax:
                    ax.scatter(punto.x_utm, punto.y_utm, marker="*", color="black",
                               edgecolor="white", s=180, linewidth=.7, zorder=11)
                    ax.annotate("KP", (punto.x_utm, punto.y_utm), xytext=(6, 8),
                                textcoords="offset points", fontsize=9)

            ax.set(title=titulo, xlim=(xmin, xmax), ylim=(ymin, ymax),
                   xlabel="Este UTM 18S (m)", ylabel="Norte UTM 18S (m)")
            ax.set_aspect("equal")
            ax.ticklabel_format(style="plain", useOffset=False)
            ax.tick_params(axis="x", rotation=30)
            ax.grid(alpha=.15)

        leyenda = [
            Line2D([0], [0], color="#124B75", lw=2, label="Río cargado desde el outlet recortado"),
            Line2D([0], [0], marker="s", color="black", ls="none", label="Outlet recortado"),
        ]
        if fin_max is not None:
            leyenda.append(Line2D([0], [0], color="#00A6B2", lw=3,
                                  label="Outlet → límite de calibración más lejano"))
        if mostrar_terrazas:
            leyenda.append(Line2D([0], [0], marker="D", color="#7943A2", ls="none",
                                  label="Ubicación real del punto de terraza"))
        if mostrar_dataciones:
            leyenda.append(Line2D([0], [0], marker="^", color="#E57A16", ls="none",
                                  label="Ubicación real de la muestra datada"))
        if not manuales.empty:
            leyenda.append(Line2D([0], [0], marker="*", color="black", ls="none",
                                  markersize=11, label="KP manual"))
        for _, limite in tabla_limites.iterrows():
            texto = (f"{limite.limite}: {limite.referencias} · cota {limite.cotas_m} m; "
                     f"cruce a {limite.distancia_m / 1000:.2f} km")
            leyenda.append(Line2D([0], [0], marker="x", color="#B52A35", ls="none",
                                  label=textwrap.fill(texto, width=48 if ancho < 10 else 70)))
        leyenda.extend(Patch(facecolor=c, edgecolor="#777777", label=et)
                       for et, c in unidades.items())
        if mostrar_litologia and rid in RIOS_SIN_CARTA_CONFIRMADA:
            leyenda.append(Patch(facecolor="white", edgecolor="#777777",
                                 label="Sin carta geológica disponible para este río"))

        ncol = 2
        # Reserva espacio para la leyenda, incluyendo etiquetas de dos líneas.
        grupos_leyenda = [leyenda[i:i + int(np.ceil(len(leyenda) / ncol))]
                         for i in range(0, len(leyenda), int(np.ceil(len(leyenda) / ncol)))]
        alto_leyenda = .45 + max(
            sum(.18 + .13 * h.get_label().count("\n") for h in grupo)
            for grupo in grupos_leyenda
        )
        alto = 6.5 + alto_leyenda
        fig.set_size_inches(ancho, alto)
        fig.subplots_adjust(left=.09 if len(cajas) == 1 else .06, right=.97,
                            top=1 - .6 / alto, bottom=(alto_leyenda + .25) / alto,
                            wspace=.3)
        subtitulo = "EPSG:32718 · Puntos de referencia reales; L = cruces de cotas sobre el río"
        if mostrar_litologia and rid not in RIOS_SIN_CARTA_CONFIRMADA:
            subtitulo += "\nBlanco: sin polígonos de la capa en ese lugar"
        fig.suptitle(f"{rid} · {NOMBRES[rid]} — {etiqueta_archivo.replace('_', ' ')}\n{subtitulo}",
                     fontsize=11)
        fig.legend(handles=leyenda, loc="lower center", bbox_to_anchor=(.5, .015),
                   fontsize=8, ncol=ncol, frameon=False)
        if guardar:
            carpeta = Path(SALIDA)
            carpeta.mkdir(parents=True, exist_ok=True)
            ruta = carpeta / f"mapa_{etiqueta_archivo}_R{rid}_{sello}.png"
            fig.savefig(ruta, dpi=300, bbox_inches="tight")
            guardadas.append(str(ruta))
        plt.show()
        plt.close(fig)
        return referencias_visibles

    if mostrar_rio_completo:
        # Incluye todo el cauce y los puntos seleccionados, aunque estén fuera del cauce.
        x = np.r_[este, puntos.este_m.to_numpy(float)]
        y = np.r_[norte, puntos.norte_m.to_numpy(float)]
        margen = max(500., .06 * max(np.ptp(x), np.ptp(y)))
        bounds = (x.min() - margen, y.min() - margen, x.max() + margen, y.max() + margen)
        dibujar_lamina([("Río completo: outlet recortado → cabecera", bounds)], "rio_completo")

    if mostrar_detalles:
        cajas = []
        if fin_max is not None:
            m = margen_calibracion_m
            cajas.append(("Detalle de calibración", (
                xcal.min() - m, ycal.min() - m, xcal.max() + m, ycal.max() + m,
            )))
        for _, punto in manuales.iterrows():
            x, y, r = float(punto.x_utm), float(punto.y_utm), radio_kp_m
            cajas.append((f"Entorno del KP · {punto.id_kp}", (x-r, y-r, x+r, y+r)))
        if cajas:
            visibles = dibujar_lamina(cajas, "detalles_calibracion_KP")
            fuera = [r.nombre for _, r in puntos.iterrows() if (r.tipo, r.nombre) not in visibles]
            if fuera:
                print("Fuera de los detalles; consulta el mapa del río completo:", ", ".join(fuera))
        else:
            print("Sin límites calculables ni KP manual para mostrar detalles.")

    print("L1, L2... identifican cruces sobre el cauce; no son las ubicaciones reales de las terrazas.")
    if not tabla_limites.empty:
        display(tabla_limites[["limite", "referencias", "cotas_m", "distancia_m"]]
                .style.format({"distancia_m": "{:.1f}"}))
    if guardar:
        carpeta = Path(SALIDA)
        carpeta.mkdir(parents=True, exist_ok=True)
        tabla_limites.to_csv(carpeta / f"limites_mapa_R{rid}_{sello}.csv",
                             index=False, encoding="utf-8-sig")
        puntos.to_csv(carpeta / f"puntos_referencia_mapa_R{rid}_{sello}.csv",
                      index=False, encoding="utf-8-sig")

    # Conservar la consulta bilateral de la litología de los KP manuales.
    filas_asignacion = []
    for _, punto in manuales.iterrows():
        asignacion = consultar_litologia(punto.x_utm, punto.y_utm, rid=rid,
                                        posicion_m=float(punto.distancia_m))
        asignacion.insert(0, "id_kp", punto.id_kp)
        asignacion.insert(0, "ID", rid)
        filas_asignacion.append(asignacion)
        print(f"KP: {punto.id_kp}")
        display(asignacion)
    if guardar and filas_asignacion:
        pd.concat(filas_asignacion, ignore_index=True).to_csv(
            Path(SALIDA) / f"litologia_KP_R{rid}_{sello}.csv",
            index=False, encoding="utf-8-sig",
        )
    if guardadas:
        print("Mapas guardados:", *guardadas, sep="\n")
    return tabla_limites


# K8 — SENSIBILIDAD A LOS CRUCES (añadir al final de K8)
# Requiere K1–K4 y las funciones originales de K8 ya definidas.
# No modifica resultados_K, CRUCES_ELEGIDOS, perfiles ni las muestras.
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

def _ajuste_intervalo_k8(p, inicio, fin):
    """Mismo OLS z–χ con intercepto libre, ahora entre dos límites."""
    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    chi = p.chi_040_m.to_numpy(float)
    if not np.isfinite([inicio, fin]).all() or not s[0] <= inicio < fin <= s[-1]:
        raise ValueError("Intervalo fuera del perfil o de longitud nula.")
    dentro = (s >= inicio) & (s <= fin)
    if int(dentro.sum()) < MIN_NODOS_K:
        raise ValueError("Pocos nodos originales para ajustar este tramo.")
    izq = max(0, int(np.searchsorted(s, inicio, side="right")) - 1)
    der = min(len(s)-1, int(np.searchsorted(s, fin, side="left")))
    if not np.isfinite(z[izq:der+1]).all() or not np.isfinite(chi[izq:der+1]).all():
        raise ValueError("El intervalo contiene datos faltantes.")
    posiciones = np.unique(np.r_[inicio, s[dentro], fin])
    x, y = np.interp(posiciones, s, chi), np.interp(posiciones, s, z)
    xc, yc = x - x.mean(), y - y.mean()
    den = float(xc @ xc)
    if den <= 0:
        raise ValueError("Sin variación de χ.")
    m = float(xc @ yc / den)
    if not np.isfinite(m) or m <= 0:
        raise ValueError("Mχ no positivo.")
    b = float(y.mean() - m*x.mean())
    residuo = y - (b + m*x)
    total = float(yc @ yc)
    return dict(Mchi=m, intercepto_m=b, ksn_040=m*A0**THETA,
                R2=1-float(residuo @ residuo)/total if total > 0 else np.nan,
                RMSE_m=float(np.sqrt(np.mean(residuo**2))))

from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import is_color_like


def _preparar_litologia_perfil_k8(p, rid):
    """Una consulta por segmento; se reutiliza para todas las terrazas del río."""
    s = p.distancia_m.to_numpy(float)
    if len(s) < 2 or not np.all(np.diff(s) > 0):
        raise ValueError("El perfil debe tener distancias crecientes y al menos dos nodos.")
    medios = (s[:-1] + s[1:]) / 2
    x = np.interp(medios, s, p.x_utm.to_numpy(float))
    y = np.interp(medios, s, p.y_utm.to_numpy(float))
    puntos = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(x, y), crs="EPSG:32718"
    )
    # Conserva intersección directa y el control de ambas riberas de K8.
    lito = asignar_unidades(puntos, ids=str(rid), posiciones_m=medios)
    lito = lito.reset_index(drop=True).copy()
    lito["s0"] = s[:-1]
    lito["s1"] = s[1:]

    colores, etiquetas = [], []
    for _, fila in lito.iterrows():
        if fila.estado_geo == "asignado":
            color = fila.color_hex
            if not isinstance(color, str) or not is_color_like(color):
                color = "#B7B7B7"
            etiqueta = f"{fila.codigo} · {fila.unidad} · {fila.fuente}"
        else:
            color = "#B7B7B7"
            if str(rid) in RIOS_SIN_CARTA_CONFIRMADA:
                etiqueta = "Sin carta geológica disponible"
            elif fila.estado_geo == "Ambiguo":
                etiqueta = "Asignación ambigua: sin unidad atribuida"
            else:
                etiqueta = "Sin unidad cartográfica asignada"
        colores.append(color)
        etiquetas.append(etiqueta)
    lito["color_grafico"] = colores
    lito["etiqueta_grafico"] = etiquetas
    return lito


def _dibujar_perfil_k8(ax, p, inicio, fin, eje="chi", lito=None):
    """Dibuja el perfil recibido. El color no interviene en los ajustes de K."""
    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    coord = p.chi_040_m.to_numpy(float) if eje == "chi" else s / 1000
    if lito is None:
        posiciones = np.unique(np.r_[inicio, s[(s > inicio) & (s < fin)], fin])
        ax.plot(np.interp(posiciones, s, coord), np.interp(posiciones, s, z),
                color="0.45", lw=1.5, label="Perfil observado", zorder=2)
        return {}

    seleccion = lito.loc[(lito.s1 > inicio) & (lito.s0 < fin)].copy()
    if seleccion.empty:
        return {}
    a = np.maximum(seleccion.s0.to_numpy(float), inicio)
    b = np.minimum(seleccion.s1.to_numpy(float), fin)
    puntos_a = np.column_stack([np.interp(a, s, coord), np.interp(a, s, z)])
    puntos_b = np.column_stack([np.interp(b, s, coord), np.interp(b, s, z)])
    segmentos = np.stack([puntos_a, puntos_b], axis=1)

    # Borde oscuro: permite ver incluso los colores claros de la carta.
    ax.add_collection(LineCollection(segmentos, colors="#333333",
                                    linewidths=4.5, zorder=1))
    ax.add_collection(LineCollection(segmentos,
                                    colors=seleccion.color_grafico.tolist(),
                                    linewidths=3.1, zorder=2))
    ax.autoscale_view()
    return dict(zip(seleccion.etiqueta_grafico, seleccion.color_grafico))


def _suavizar_perfil_k8(p, ventana_m=100.0, metodo="savitzky"):
    """Suaviza la elevación sobre una malla regular en distancia.

    metodo="savitzky": ajuste local de grado 2.
    metodo="media_movil": promedio aritmético centrado, pesos iguales.
    Ambos usan la misma regla de ventana impar y redondeo hacia arriba.
    """
    from scipy.signal import savgol_filter
    from scipy.ndimage import uniform_filter1d

    metodo = str(metodo).strip().lower()
    if metodo not in {"savitzky", "media_movil"}:
        raise ValueError('Método no válido: elige "savitzky" o "media_movil".')

    s = p.distancia_m.to_numpy(float)
    z = p.elevacion_m.to_numpy(float)
    if len(s) < 5 or not np.isfinite(s).all() or not np.isfinite(z).all():
        raise ValueError("El suavizado requiere al menos cinco nodos sin datos faltantes.")
    if not np.all(np.diff(s) > 0):
        raise ValueError("Las distancias del perfil deben ser crecientes.")
    ventana_m = float(ventana_m)
    if not np.isfinite(ventana_m) or ventana_m <= 0:
        raise ValueError("La ventana de suavizado debe ser positiva y finita.")

    paso_aprox = float(np.median(np.diff(s)))
    n = int(np.ceil((s[-1] - s[0]) / paso_aprox)) + 1
    malla = np.linspace(s[0], s[-1], n)
    paso = float(malla[1] - malla[0])

    # Número impar de nodos; registra la anchura realmente utilizada.
    nodos_ventana = int(np.ceil(ventana_m / paso)) + 1
    if nodos_ventana % 2 == 0:
        nodos_ventana += 1
    if nodos_ventana < 5:
        raise ValueError(
            f"Ventana demasiado pequeña: usa al menos {4 * paso:.1f} m."
        )
    if nodos_ventana > n:
        raise ValueError("La ventana de suavizado supera la longitud del perfil.")

    z_regular = np.interp(malla, s, z)
    if metodo == "savitzky":
        z_filtrado = savgol_filter(
            z_regular, window_length=nodos_ventana,
            polyorder=2, mode="interp",
        )
        nombre_metodo = "Savitzky–Golay"
        modo_borde = "Ajuste polinómico en extremos"
    else:
        # Promedio aritmético de los nodos de una ventana centrada.
        # En los extremos, refleja los datos para completar la ventana.
        z_filtrado = uniform_filter1d(
            z_regular, size=nodos_ventana, mode="reflect", origin=0,
        )
        nombre_metodo = "Media móvil"
        modo_borde = "Reflexión en extremos"
    p_suavizado = p.copy(deep=True)
    p_suavizado["elevacion_m"] = np.interp(s, malla, z_filtrado)

    # Conserva distancia, chi, area y coordenadas originales.
    return p_suavizado, {
        "ventana_solicitada_m": ventana_m,
        "ventana_real_m": (nodos_ventana - 1) * paso,
        "paso_malla_m": paso,
        "nodos_ventana": nodos_ventana,
        "orden_polynomial": 2 if metodo == "savitzky" else np.nan,
        "metodo_suavizado": nombre_metodo,
        "metodo_clave": metodo,
        "tratamiento_extremos": modo_borde,
    }


def sensibilidad_cruces_k8(rid="9", referencias=None, muestras=None,
                           id_kp=None, mostrar_mapa=False, guardar=False,
                           mostrar_litologia=True, suavizado=False,
                           ventana_suavizado_m=100.0,
                           mapa_rio_completo=True, mapa_detalles=True,
                           mapa_terrazas=True, mapa_dataciones=True, ver_juntos=False,
                           metodo_suavizado="savitzky"):
    """Cada modelo usa sus elevaciones, sus cruces y sus ajustes de K.

    Raw: original. Suavizado: método elegido sobre elevación–distancia.
    Se conservan χ, área, U, cota de referencia y posición espacial del KP.
    El primer cruce desde el outlet es principal dentro de cada modelo.
    No modifica perfiles, resultados_K ni CRUCES_ELEGIDOS.
    """
    requeridas = ["calibracion", "perfiles", "buscar_cruces", "ajustar_hasta_cruce",
                  "kp_manuales", "MIN_NODOS_K", "THETA", "A0", "N", "RIOS_K", "NOMBRES"]
    faltan = [nombre for nombre in requeridas if nombre not in globals()]
    if faltan:
        raise RuntimeError(f"Ejecuta las celdas anteriores. Faltan: {faltan}")
    if not (np.isclose(THETA, .40) and np.isclose(N, 1)):
        raise ValueError("Se requiere θ=0,40 y n=1.")
    rid = str(rid)
    if rid not in perfiles or rid not in list(map(str, RIOS_K)):
        raise ValueError("Usa uno de los seis ríos cargados.")
    seleccion = calibracion.loc[
        calibracion.ID.astype(str).eq(rid) & calibracion.tipo_cota.eq("terraza")
    ].copy()
    for columna, filtro in [("referencia", referencias), ("muestra", muestras)]:
        if filtro is not None:
            filtro = [filtro] if isinstance(filtro, str) else list(filtro)
            if filtro:
                desconocidos = set(filtro) - set(seleccion[columna])
                if desconocidos:
                    raise ValueError(f"{columna} no disponible: {desconocidos}")
                seleccion = seleccion.loc[seleccion[columna].isin(filtro)]
    if seleccion.empty:
        raise ValueError("No hay escenarios para esa selección.")

    puntos_kp = kp_manuales.loc[kp_manuales.ID.astype(str).eq(rid)].copy()
    if id_kp is not None:
        puntos_kp = puntos_kp.loc[puntos_kp.id_kp.eq(id_kp)]
        if puntos_kp.empty:
            raise ValueError("El KP indicado no existe en la selección manual de este río.")
    if len(puntos_kp) > 1:
        raise ValueError("Hay varios KP manuales: indica id_kp explícitamente.")
    punto_kp = None if puntos_kp.empty else puntos_kp.iloc[0]

    p = perfiles[rid]
    s = p.distancia_m.to_numpy(float)
    chi = p.chi_040_m.to_numpy(float)
    fin_kp = None if punto_kp is None else float(punto_kp.distancia_m)
    lito_perfil = _preparar_litologia_perfil_k8(p, rid) if mostrar_litologia else None
    modo = "ambos" if ver_juntos else ("suavizado" if suavizado else "raw")
    modelos = {}
    if ver_juntos or not suavizado:
        modelos["raw"] = {"perfil": p, "suavizado": "No", "info": {}}
    if ver_juntos or suavizado:
        p_suavizado, info = _suavizar_perfil_k8(p, ventana_suavizado_m, metodo=metodo_suavizado)
        modelos["suavizado"] = {"perfil": p_suavizado, "suavizado": "Sí", "info": info}
        print(f"{info['metodo_suavizado']}: ventana solicitada {info['ventana_solicitada_m']:.1f} m; "
              f"ventana efectiva {info['ventana_real_m']:.1f} m.")
    print(f"Vista: {modo}. Cada curva usa sus propios cruces y ajustes.")

    def intentar_ajuste(perfil, inicio, fin, es_k1=False):
        try:
            ajuste = (ajustar_hasta_cruce(perfil, fin) if es_k1
                      else _ajuste_intervalo_k8(perfil, inicio, fin))
            return ajuste, "Calculado condicional"
        except ValueError as exc:
            return None, str(exc)

    todas = []
    sello = datetime.now(ZoneInfo("America/Los_Angeles")).strftime("%Y%m%d_%H%M%S_%f")
    for _, escenario in seleccion.iterrows():
        u = float(escenario.U_m_ka)
        cota = float(escenario.cota_m)
        if not np.isfinite([u, cota]).all() or u <= 0:
            raise ValueError("U debe ser positiva y finita, y la cota debe ser finita.")
        datos_modelos, filas = {}, []
        for nombre_modelo, modelo in modelos.items():
            perfil_modelo = modelo["perfil"]
            # La búsqueda se repite sobre la curva elegida: no se heredan
            # los cruces del original para el modelo suavizado.
            cruces, planos = buscar_cruces(perfil_modelo, cota)
            cruces = sorted(float(x) for x in cruces)
            ajuste3, estado3 = None, "Sin KP manual"
            if fin_kp is not None:
                ajuste3, estado3 = intentar_ajuste(perfil_modelo, s[0], fin_kp)
            datos_modelos[nombre_modelo] = {
                **modelo, "cruces": cruces, "ajustes": {}, "ajuste3": ajuste3,
            }

            # Una fila por cruce y modelo. Si no hay cruce, conservar
            # el escenario como pendiente; K3 puede existir igualmente.
            for numero, fin in enumerate(cruces if cruces else [np.nan], 1):
                fila = dict(
                    escenario_id=escenario.escenario_id,
                    ID=rid, rio=NOMBRES[rid], referencia=escenario.referencia,
                    muestra=escenario.muestra, cota_m=cota, U_m_ka=u,
                    modelo=nombre_modelo, suavizado=modelo["suavizado"],
                    cruce=numero if cruces else np.nan, n_cruces=len(cruces),
                    principal=bool(cruces and numero == 1), distancia_m=fin,
                    id_kp=None if punto_kp is None else punto_kp.id_kp,
                    origen_cruces="Curva suavizada" if nombre_modelo == "suavizado" else "DEM original",
                )
                if nombre_modelo == "suavizado":
                    fila.update(modelo["info"])
                else:
                    fila.update(metodo_suavizado="Sin suavizado", metodo_clave="raw",
                                ventana_solicitada_m=0.0, ventana_real_m=0.0,
                                tratamiento_extremos="Sin filtro")
                a1, a2 = None, None
                estado1, estado2 = "Sin cruces", "Sin KP manual" if fin_kp is None else "Sin cruces"
                if cruces:
                    plano = any(a - 1e-6 <= fin <= b + 1e-6 for a, b in planos)
                    if plano:
                        estado1 = "Cruce plano: revisar"
                        if fin_kp is not None:
                            estado2 = "Cruce plano: revisar"
                    else:
                        a1, estado1 = intentar_ajuste(perfil_modelo, s[0], fin, es_k1=True)
                        if fin_kp is not None:
                            a2, estado2 = intentar_ajuste(perfil_modelo, fin, fin_kp)
                    datos_modelos[nombre_modelo]["ajustes"][numero] = (a1, a2)
                for campo, ajuste, estado in [
                    ("K1", a1, estado1), ("K2", a2, estado2), ("K3", ajuste3, estado3),
                ]:
                    fila[campo] = np.nan if ajuste is None else (u / 1000) / ajuste["ksn_040"]
                    fila[f"estado_{campo}"] = estado
                    for metrica in ["Mchi", "R2", "RMSE_m"]:
                        fila[f"{metrica}_{campo}"] = np.nan if ajuste is None else ajuste[metrica]
                filas.append(fila)

        tabla = pd.DataFrame(filas)
        todas.extend(filas)
        print(f"\n{NOMBRES[rid]} · {escenario.referencia} · {cota:.2f} m · {escenario.muestra}")
        print("Principal: primer cruce desde el outlet, por separado en cada modelo.")
        columnas = ["suavizado", "metodo_suavizado", "ventana_real_m",
                    "cruce", "principal", "distancia_m", "K1", "K2", "K3"]
        display(tabla[columnas].rename(columns={
            "suavizado": "Suavizado", "metodo_suavizado": "Método",
            "ventana_real_m": "Ventana efectiva (m)", "cruce": "Cruce", "principal": "Principal",
            "distancia_m": "Distancia del cruce (m)",
        }).style.format({
            "Cruce": "{:.0f}", "Distancia del cruce (m)": "{:.1f}",
            "Ventana efectiva (m)": "{:.1f}",
            "K1": "{:.3e}", "K2": "{:.3e}", "K3": "{:.3e}",
        }, na_rep="—"))
        estados = ["estado_K1", "estado_K2", "estado_K3"]
        pendientes = tabla.loc[
            ~tabla[estados].isin(["Calculado condicional", "Sin KP manual"]).all(axis=1)
        ]
        if not pendientes.empty:
            display(pendientes[["suavizado", "cruce", *estados]])

        fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
        leyenda_lito = {}
        todos_cruces = [x for datos in datos_modelos.values() for x in datos["cruces"]]
        if todos_cruces:
            margen = max(100., .2 * (max(todos_cruces) - min(todos_cruces)))
            lo, hi = max(s[0], min(todos_cruces) - margen), min(s[-1], max(todos_cruces) + margen)
        else:
            lo, hi = s[0], s[-1]
        extremo = max([s[0], *todos_cruces, *([] if fin_kp is None else [fin_kp])])
        hasta = min(s[-1], extremo + max(100., .04 * extremo)) if extremo > s[0] else s[-1]

        for nombre_modelo, datos in datos_modelos.items():
            es_suavizado = nombre_modelo == "suavizado"
            perfil_modelo = datos["perfil"]
            z_modelo = perfil_modelo.elevacion_m.to_numpy(float)
            color_modelo = "black" if es_suavizado else "0.55"
            etiqueta_modelo = "Suavizado" if es_suavizado else "Raw"
            marcador = "s" if es_suavizado else "o"
            estilo = ":" if es_suavizado and ver_juntos else "--"
            for eje_ax, inicio, fin, eje in [
                (ax[0], lo, hi, "distancia"), (ax[1], s[0], hasta, "chi"),
            ]:
                if ver_juntos and not es_suavizado:
                    pos = np.unique(np.r_[inicio, s[(s > inicio) & (s < fin)], fin])
                    xx = pos / 1000 if eje == "distancia" else np.interp(pos, s, chi)
                    eje_ax.plot(xx, np.interp(pos, s, z_modelo), color=color_modelo,
                                lw=1.3, label="Perfil raw", zorder=3)
                else:
                    antes = len(eje_ax.lines)
                    leyenda_lito.update(_dibujar_perfil_k8(
                        eje_ax, perfil_modelo, inicio, fin, eje=eje, lito=lito_perfil,
                    ))
                    if len(eje_ax.lines) > antes:
                        eje_ax.lines[-1].set_label(f"Perfil {etiqueta_modelo.lower()}")
                        eje_ax.lines[-1].set_color(color_modelo)
                    else:
                        eje_ax.plot([], [], color=color_modelo, label=f"Perfil {etiqueta_modelo.lower()}")
            for numero, cruce in enumerate(datos["cruces"], 1):
                color = "#1261a0" if numero == 1 else "#e28b24"
                for eje_ax, xx in [(ax[0], cruce / 1000), (ax[1], np.interp(cruce, s, chi))]:
                    eje_ax.scatter(xx, cota, marker=marcador, color=color, s=50,
                                   edgecolors="white", linewidths=.5, zorder=6,
                                   label=f"{etiqueta_modelo} · {'primer cruce' if numero == 1 else 'otros cruces'}"
                                   if numero <= 2 else None)
                ax[0].annotate(f"{etiqueta_modelo} {numero}", (cruce / 1000, cota),
                               xytext=(0, (10 + 12 * ((numero - 1) % 3)) * (-1 if es_suavizado else 1)),
                               textcoords="offset points", ha="center", fontsize=7)

            if datos["cruces"]:
                fin = datos["cruces"][0]
                a1, a2 = datos["ajustes"][1]
                for ajuste, inicio, final, color, etiqueta in [
                    (a1, s[0], fin, "#1261a0", "K1"),
                    (a2, fin, fin_kp, "#d46c00", "K2"),
                    (datos["ajuste3"], s[0], fin_kp, "0.25", "K3"),
                ]:
                    if ajuste is not None:
                        xx = np.interp([inicio, final], s, chi)
                        ax[1].plot(xx, ajuste["intercepto_m"] + ajuste["Mchi"] * xx,
                                   color=color, ls=estilo, lw=1.8, zorder=4,
                                   label=f"{etiqueta} {etiqueta_modelo}: Mχ={ajuste['Mchi']:.4f}")
            elif datos["ajuste3"] is not None:
                ajuste = datos["ajuste3"]
                xx = np.interp([s[0], fin_kp], s, chi)
                ax[1].plot(xx, ajuste["intercepto_m"] + ajuste["Mchi"] * xx,
                           color="0.25", ls=estilo, label=f"K3 {etiqueta_modelo}: Mχ={ajuste['Mchi']:.4f}")
            if fin_kp is not None:
                # Misma ubicación manual; elevación sobre la curva de cada modelo.
                ax[1].scatter(np.interp(fin_kp, s, chi), np.interp(fin_kp, s, z_modelo),
                              marker="*", color=color_modelo, edgecolors="white", s=140,
                              linewidths=.6, zorder=7, label=f"KP manual · {etiqueta_modelo}")
            subtabla = tabla.loc[tabla.modelo.eq(nombre_modelo) & tabla.cruce.notna()]
            for campo, color in [("K1", "#1261a0"), ("K2", "#d46c00"), ("K3", "0.3")]:
                if subtabla[campo].notna().any():
                    ax[2].plot(subtabla.cruce, subtabla[campo], marker=marcador,
                               ls=estilo, color=color, label=f"{campo} {etiqueta_modelo}")

        ax[0].axhline(cota, color="#9b4a00", ls="--",
                      label=f"{escenario.referencia} · {cota:.2f} m")
        ax[0].set(title="Cruces propios de cada curva", xlabel="Distancia desde outlet (km)",
                  ylabel="Elevación (m)", xlim=(lo / 1000, hi / 1000))
        ax[1].set(title="Ajustes del primer cruce de cada modelo", xlabel="χ (m)", ylabel="Elevación (m)")
        ax[2].set(title="Sensibilidad al cruce, por modelo", xlabel="Número de cruce dentro de cada modelo",
                  ylabel="K (m⁰·²/año)")
        ax[2].set_xticks(sorted(tabla.cruce.dropna().unique()))
        ax[2].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        for eje_ax in ax:
            eje_ax.grid(alpha=.2)
            handles, labels = eje_ax.get_legend_handles_labels()
            if handles:
                leyenda = dict(zip(labels, handles))
                eje_ax.legend(leyenda.values(), leyenda.keys(), fontsize=7)
        detalle = "Raw"
        if "suavizado" in modelos:
            detalle = ("Raw + suavizado" if ver_juntos else "Suavizado")
            detalle += (f" · {modelos['suavizado']['info']['metodo_suavizado']}"
                        f" {modelos['suavizado']['info']['ventana_real_m']:.1f} m")
        fig.suptitle(f"{rid} · {NOMBRES[rid]} · {escenario.referencia} · {escenario.muestra}\n"
                     f"{detalle} · Cruces recalculados por modelo", fontsize=11)
        if mostrar_litologia and leyenda_lito:
            elementos = [Patch(facecolor=color, edgecolor="#333333", label=etiqueta)
                         for etiqueta, color in leyenda_lito.items()]
            espacio = .55 + .24 * len(elementos)
            altura = 4.8 + espacio
            fig.set_size_inches(17, altura)
            fig.subplots_adjust(left=.06, right=.98, top=1 - .55 / altura,
                                bottom=espacio / altura, wspace=.35)
            fig.legend(handles=elementos, title="Litología cartográfica", loc="lower center",
                       bbox_to_anchor=(.5, .015), ncol=1, fontsize=8, title_fontsize=9, frameon=False)
        else:
            fig.tight_layout(rect=(0, 0, 1, .93))
        if guardar:
            carpeta = Path(SALIDA)
            carpeta.mkdir(parents=True, exist_ok=True)
            sg = (f"_{modelos['suavizado']['info']['metodo_clave']}"
                  f"_{modelos['suavizado']['info']['ventana_real_m']:.1f}m"
                  if "suavizado" in modelos else "")
            lito = "con_litologia" if mostrar_litologia else "sin_litologia"
            etiqueta = f"R{rid}_{escenario.referencia}_{escenario.muestra.replace(' ', '')}_{modo}{sg}_{sello}"
            fig.savefig(carpeta / f"cruces_modelos_{lito}_{etiqueta}.png", dpi=300, bbox_inches="tight")
        plt.show()
        plt.close(fig)

    resultado = pd.DataFrame(todas)
    if mostrar_mapa:
        # El mapa usa ahora los primeros cruces de los modelos visibles.
        limites_mapa = resultado.loc[resultado.principal & resultado.distancia_m.notna()].copy()
        limites_mapa["fin_m"] = limites_mapa.distancia_m
        limites_mapa["referencia"] = limites_mapa.referencia.astype(str) + " · " + limites_mapa.modelo
        mapa_revision(rid, guardar=guardar, mostrar_rio_completo=mapa_rio_completo,
                      mostrar_detalles=mapa_detalles, mostrar_terrazas=mapa_terrazas,
                      mostrar_dataciones=mapa_dataciones, mostrar_litologia=mostrar_litologia,
                      referencias=referencias, muestras=muestras, limites_modelos=limites_mapa)
    if guardar:
        carpeta = Path(SALIDA)
        carpeta.mkdir(parents=True, exist_ok=True)
        resultado.to_csv(carpeta / f"K_cruces_modelos_R{rid}_{modo}_{sello}.csv",
                         index=False, encoding="utf-8-sig", float_format="%.12g")
        print("Tabla y figuras guardadas en:", carpeta)
    print("K1: outlet–cruce de cota; K2: cruce–KP; K3: outlet–KP.")
    print("K3 se calcula una vez por modelo y escenario; se repite para acompañar cada cruce.")
    return resultado


# K8 — CONTROLES Y EJECUCIÓN
RIOS_VER_K8 = ["10"]  # [] incluye los seis.
REFERENCIAS_VER_K8 = []  # [] incluye todas las terrazas de cada río.
MUESTRAS_VER_K8 = []     # [] incluye todas las muestras de cada río.

MOSTRAR_LITOLOGIA_K8 = True  # False dibuja los perfiles en gris.
SUAVIZADO_K8 = False         # True: solo suavizado. False: solo original.
VER_JUNTOS_K8 = False       # True: raw + el método de suavizado elegido.

# -------- ELEGIR EL MÉTODO DE SUAVIZADO --------
METODO_SUAVIZADO_K8 = "savitzky"   # Opciones: "savitzky" o "media_movil"
VENTANA_SUAVIZADO_K8_M = 100.0  # Metros a lo largo del cauce.
# Activar el mapa y elegir sus vistas/contenido:
MOSTRAR_MAPA_K8 = True
MOSTRAR_RIO_COMPLETO_K8 = True
MOSTRAR_DETALLES_MAPA_K8 = True
MOSTRAR_TERRAZAS_MAPA_K8 = True
MOSTRAR_DATACIONES_MAPA_K8 = True

# Coordenadas verificadas del Ages_Final actual (encabezados invertidos).
COLUMNA_ESTE_DATACION_K8 = "Ny"
COLUMNA_NORTE_DATACION_K8 = "Ex"

GUARDAR_K8 = False          # True guarda tablas/figuras con fecha y hora.

rios_k8 = [str(rid) for rid in RIOS_VER_K8] if RIOS_VER_K8 else ["3", "4", "7", "8", "9", "10"]
tablas_k8 = []
for rid in rios_k8:
    resultado_rio = sensibilidad_cruces_k8(
        rid=rid,
        referencias=REFERENCIAS_VER_K8,
        muestras=MUESTRAS_VER_K8,
        mostrar_mapa=MOSTRAR_MAPA_K8,
        guardar=GUARDAR_K8,
        mostrar_litologia=MOSTRAR_LITOLOGIA_K8,
        suavizado=SUAVIZADO_K8,
        ver_juntos=VER_JUNTOS_K8,
        ventana_suavizado_m=VENTANA_SUAVIZADO_K8_M,
        metodo_suavizado=METODO_SUAVIZADO_K8,
        mapa_rio_completo=MOSTRAR_RIO_COMPLETO_K8,
        mapa_detalles=MOSTRAR_DETALLES_MAPA_K8,
        mapa_terrazas=MOSTRAR_TERRAZAS_MAPA_K8,
        mapa_dataciones=MOSTRAR_DATACIONES_MAPA_K8,
    )
    tablas_k8.append(resultado_rio)

tabla_sensibilidad_K8 = pd.concat(tablas_k8, ignore_index=True)


# Notebook cell 22

