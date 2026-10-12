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

pa.NOTEBOOK_PATH = CARPETA_NOTEBOOKS / 'AnalisisPerfiles_20261005.ipynb'


# Notebook cell 3
from scipy.integrate import cumulative_trapezoid
from scipy.stats import linregress
from scipy.signal import find_peaks
import hashlib
perfiles_originales, nombres, manifest = pa.load_profiles()
deps = pa.dependency('perfiles')
deps[str(Path(CFG['kp_manual_csv']))] = pa.file_hash(CFG['kp_manual_csv'])
pa.begin('analisis_perfiles')
BASE = Path(CFG['base_dir']) / 'Productos'
EJECUTAR_LOCAL = True
EJECUTAR_TOPOTOOLBOX = False  # Requiere instalar/configurar el entorno ya usado por el original.
EJECUTAR_SENSIBILIDAD_LOCAL = False
EJECUTAR_SENSIBILIDAD_TT = False


# Notebook cell 5
# ── MODIFICAR AQUÍ: intervalos que quieres analizar ──────────────────────
# Formato: "ID": [("nombre del tramo", inicio_km, fin_km), ...]

TRAMOS = {
    # "6": [
    #     ("inferior", 2.0, 8.0),
    #     ("superior", 8.0, None),
    # ],
}

# Los intervalos son [inicio, fin): el extremo final no se duplica.
# Para excluir una zona intermedia, deja un espacio entre intervalos.

ids_desconocidos = set(TRAMOS) - set(perfiles_originales)
if ids_desconocidos:
    raise ValueError(f"IDs desconocidos en TRAMOS: {ids_desconocidos}")

perfiles = {}
mascaras_tramos = {}
registros_qc = []

for river_id, original in perfiles_originales.items():
    p = original.copy(deep=True)

    x = p["distancia_m"].to_numpy(dtype=float)
    z = p["elevacion_m"].to_numpy(dtype=float)
    a = p["area_m2"].to_numpy(dtype=float)

    # No integrar sobre distancias o áreas inválidas
    if len(p) < 2 or not np.all(np.isfinite(x)):
        raise ValueError(f"{river_id}: distancia insuficiente o inválida.")

    if np.any(np.diff(x) <= 0):
        raise ValueError(f"{river_id}: la distancia no aumenta estrictamente.")

    if not np.isclose(x[0], 0):
        raise ValueError(f"{river_id}: el perfil no comienza en el outlet.")

    if np.any(~np.isfinite(a) | (a <= 0)):
        raise ValueError(
            f"{river_id}: hay áreas inválidas. Revisar sin borrar las filas."
        )

    p["elevacion_valida"] = np.isfinite(z)

    # Registrar inversiones o zonas planas, sin corregirlas ni usar abs(S)
    pendiente_entre_puntos = np.concatenate((
        [np.nan], np.diff(z) / np.diff(x)
    ))
    p["pendiente_paso"] = pendiente_entre_puntos
    p["paso_no_positivo"] = (
        np.isfinite(pendiente_entre_puntos)
        & (pendiente_entre_puntos <= 0)
    )

    intervalos = TRAMOS.get(
        river_id, [("completo_provisional", 0.0, None)]
    )

    mascaras_tramos[river_id] = {}
    asignado = np.zeros(len(p), dtype=bool)

    for nombre_tramo, inicio_km, fin_km in intervalos:
        if nombre_tramo in mascaras_tramos[river_id]:
            raise ValueError(f"{river_id}: nombre de tramo repetido.")

        inicio = float(inicio_km) * 1000
        fin = np.inf if fin_km is None else float(fin_km) * 1000

        if inicio < 0 or fin <= inicio:
            raise ValueError(f"{river_id}: límites inválidos.")

        mask = (x >= inicio) & (x < fin)

        if not mask.any():
            raise ValueError(f"{river_id} — {nombre_tramo}: tramo vacío.")

        if np.any(mask & asignado):
            raise ValueError(f"{river_id}: hay tramos superpuestos.")

        mascaras_tramos[river_id][nombre_tramo] = mask
        asignado |= mask

    p["incluido_en_tramos"] = asignado
    perfiles[river_id] = p

    registros_qc.append({
        "ID": river_id,
        "rio": nombres[river_id],
        "longitud_km": x[-1] / 1000,
        "puntos": len(p),
        "sin_elevacion": int((~np.isfinite(z)).sum()),
        "pasos_no_positivos": int(p["paso_no_positivo"].sum()),
        "fuera_de_tramos": int((~asignado).sum()),
        "delimitacion": (
            "manual" if river_id in TRAMOS else "completa_provisional"
        ),
    })

control_calidad = pd.DataFrame(registros_qc)

display(control_calidad)


# Notebook cell 7
# ── Parámetros del análisis ─────────────────────────────────────────────
THETAS = np.arange(30, 71, 5) / 100   # Incluye 0.30 y 0.70
THETA_AJUSTE = CFG['theta']                  # Referencia inicial, ajustable
A0_ANALISIS = CFG['A0_m2']                    # m²

# Criterios iniciales para evitar ajustes demasiado cortos
MIN_PUNTOS_AJUSTE = 10
MIN_LONGITUD_M = 500.0

theta_calculados = sorted(set([
    *np.round(THETAS, 8),
    round(THETA_AJUSTE, 8),
]))

chi_sensibilidad = {}
filas_sensibilidad = []
comparacion_chi = []


def ajustar_recta(chi, z, distancia, mask):
    """Ajustar sin eliminar ni modificar filas de la tabla original."""
    indices = np.flatnonzero(mask)
    if len(indices) and (np.any(np.diff(indices) > 1) or not np.isfinite(z[indices]).all()):
        return {'estado': 'hueco_en_tramo', 'n': int(np.isfinite(z[indices]).sum())}
    valid = mask & np.isfinite(chi) & np.isfinite(z)
    n = int(valid.sum())

    if n < MIN_PUNTOS_AJUSTE:
        return {"estado": "pocos_puntos", "n": n}

    extension = float(np.ptp(distancia[valid]))
    if extension < MIN_LONGITUD_M:
        return {"estado": "tramo_corto", "n": n}

    if np.ptp(chi[valid]) == 0 or np.ptp(z[valid]) == 0:
        return {"estado": "sin_variacion", "n": n}

    fit = linregress(chi[valid], z[valid])
    residuos = z[valid] - (fit.intercept + fit.slope * chi[valid])

    return {
        "estado": "ok" if fit.slope > 0 else "pendiente_no_positiva",
        "n": n,
        "M_chi": float(fit.slope),
        "intercepto_m": float(fit.intercept),
        "R2": float(fit.rvalue**2),
        "RMSE_m": float(np.sqrt(np.mean(residuos**2))),
        "chi_inicio": float(chi[valid].min()),
        "chi_fin": float(chi[valid].max()),
    }


for river_id, p in perfiles.items():
    x = p["distancia_m"].to_numpy(dtype=float)
    a = p["area_m2"].to_numpy(dtype=float)
    z = p["elevacion_m"].to_numpy(dtype=float)

    chi_sensibilidad[river_id] = {}

    for theta in theta_calculados:
        chi = cumulative_trapezoid(
            (A0_ANALISIS / a)**theta,
            x=x,
            initial=0.0,
        )
        chi_sensibilidad[river_id][theta] = chi

        for tramo, mask in mascaras_tramos[river_id].items():
            ajuste = ajustar_recta(chi, z, x, mask)

            filas_sensibilidad.append({
                "ID": river_id,
                "rio": nombres[river_id],
                "tramo": tramo,
                "theta": theta,
                **ajuste,
            })

    # Comparación con χ previo usando SUS parámetros de referencia
    anteriores = manifest["rios"][river_id]["parametros"]

    if (
        "chi_m" in p
        and "theta_ref" in anteriores
        and "A0_m2" in anteriores
    ):
        chi_control = cumulative_trapezoid(
            (anteriores["A0_m2"] / a)**anteriores["theta_ref"],
            x=x,
            initial=0.0,
        )
        chi_previo = p["chi_m"].to_numpy(dtype=float)
        valid = np.isfinite(chi_previo)

        comparacion_chi.append({
            "ID": river_id,
            "diferencia_max_chi_m": (
                float(np.max(np.abs(
                    chi_control[valid] - chi_previo[valid]
                ))) if valid.any() else np.nan
            ),
        })

tabla_sensibilidad = pd.DataFrame(filas_sensibilidad)
#display(tabla_sensibilidad)

print("Comparación de discretizaciones de χ; no sustituye una validación:")
display(pd.DataFrame(comparacion_chi))


# Notebook cell 8
ajustes_segmentos = tabla_sensibilidad.loc[np.isclose(tabla_sensibilidad['theta'], THETA_AJUSTE)].copy()
ajustes_segmentos['A0_m2'] = A0_ANALISIS
ajustes_segmentos['ksn_segmento'] = np.where(ajustes_segmentos.estado.eq('ok'),
    ajustes_segmentos.get('M_chi', np.nan)*A0_ANALISIS**THETA_AJUSTE, np.nan)


# Notebook cell 10
from scipy.signal import find_peaks

# ── Parámetros iniciales: ajustables, no calibrados para tus ríos ─────────
THETA_DETECCION = round(float(THETA_AJUSTE), 8)

VENTANA_LADO_M = 500.0      # Longitud de cada tramo vecino
PASO_BUSQUEDA_M = 100.0     # Separación entre posiciones evaluadas
SEPARACION_CANDIDATOS_M = 500.0

MIN_PUNTOS_LADO = 10
R2_MINIMO = 0.80

# Cambio relativo = |ksn_arriba - ksn_abajo| / max(ksn_arriba, ksn_abajo)
CAMBIO_MINIMO = 0.30        # Diferencia de al menos 30 %
PROMINENCIA_MINIMA = 0.10   # Destacar frente a cambios vecinos

COLUMNAS_DETECCION = [
    "ID", "rio", "tramo", "bloque", "indice",
    "distancia_m", "x_utm", "y_utm", "elevacion_m", "chi_m",
    "theta", "ventana_lado_m",
    "ksn_abajo", "ksn_arriba", "cambio_relativo",
    "r2_abajo", "r2_arriba", "evaluacion_valida", "sentido",
]


def detectar_cambios(river_id, theta, ventana_m):
    """Detector local exploratorio; no modifica los perfiles originales."""
    theta = round(float(theta), 8)

    if theta not in chi_sensibilidad[river_id]:
        raise ValueError(f"Falta calcular θ={theta} para el río {river_id}.")

    p = perfiles[river_id]
    x = p["distancia_m"].to_numpy(dtype=float)
    z = p["elevacion_m"].to_numpy(dtype=float)
    chi = np.asarray(chi_sensibilidad[river_id][theta])

    evaluaciones = []
    candidatos = []

    for tramo, mask in mascaras_tramos[river_id].items():
        indices = np.flatnonzero(
            mask & np.isfinite(z) & np.isfinite(chi)
        )

        # No unir puntos a través de huecos sin datos
        bloques = np.split(
            indices, np.flatnonzero(np.diff(indices) > 1) + 1
        )

        for bloque_id, idx in enumerate(bloques):
            if len(idx) < 2 * MIN_PUNTOS_LADO:
                continue

            xb, zb, cb = x[idx], z[idx], chi[idx]

            if xb[-1] - xb[0] <= 2 * ventana_m:
                continue

            centros = np.arange(
                xb[0] + ventana_m,
                xb[-1] - ventana_m,
                PASO_BUSQUEDA_M,
            )

            filas_bloque = []

            for centro in centros:
                abajo = (xb >= centro - ventana_m) & (xb < centro)
                arriba = (xb > centro) & (xb <= centro + ventana_m)

                j = int(idx[np.argmin(np.abs(xb - centro))])

                fila = {
                    "ID": river_id,
                    "rio": nombres[river_id],
                    "tramo": tramo,
                    "bloque": bloque_id,
                    "indice": j,
                    "distancia_m": float(x[j]),
                    "x_utm": float(p.iloc[j]["x_utm"]),
                    "y_utm": float(p.iloc[j]["y_utm"]),
                    "elevacion_m": float(z[j]),
                    "chi_m": float(chi[j]),
                    "theta": theta,
                    "ventana_lado_m": ventana_m,
                    "ksn_abajo": np.nan,
                    "ksn_arriba": np.nan,
                    "cambio_relativo": np.nan,
                    "r2_abajo": np.nan,
                    "r2_arriba": np.nan,
                    "evaluacion_valida": False,
                    "sentido": "",
                }

                suficientes = (
                    abajo.sum() >= MIN_PUNTOS_LADO
                    and arriba.sum() >= MIN_PUNTOS_LADO
                )

                if suficientes:
                    fit_a = linregress(cb[abajo], zb[abajo])
                    fit_b = linregress(cb[arriba], zb[arriba])

                    ka = float(fit_a.slope * A0_ANALISIS**theta)
                    kb = float(fit_b.slope * A0_ANALISIS**theta)
                    r2a = float(fit_a.rvalue**2)
                    r2b = float(fit_b.rvalue**2)

                    valido = (
                        np.isfinite([ka, kb, r2a, r2b]).all()
                        and ka > 0 and kb > 0
                        and min(r2a, r2b) >= R2_MINIMO
                    )

                    fila.update({
                        "ksn_abajo": ka,
                        "ksn_arriba": kb,
                        "r2_abajo": r2a,
                        "r2_arriba": r2b,
                        "evaluacion_valida": valido,
                    })

                    if valido:
                        fila["cambio_relativo"] = abs(kb - ka) / max(ka, kb)
                        fila["sentido"] = (
                            "aumenta_aguas_arriba"
                            if kb > ka else "disminuye_aguas_arriba"
                        )

                filas_bloque.append(fila)

            if not filas_bloque:
                continue

            tabla = pd.DataFrame(filas_bloque)
            señal = np.where(
                tabla["evaluacion_valida"],
                tabla["cambio_relativo"],
                0.0,
            )

            picos, _ = find_peaks(
                señal,
                height=CAMBIO_MINIMO,
                prominence=PROMINENCIA_MINIMA,
                distance=max(
                    1,
                    int(np.ceil(
                        SEPARACION_CANDIDATOS_M / PASO_BUSQUEDA_M
                    )),
                ),
            )

            for k in picos:
                # Evitar candidatos inmediatamente junto a ajustes inválidos
                if tabla.iloc[k-1:k+2]["evaluacion_valida"].all():
                    candidatos.append(filas_bloque[k].copy())

            evaluaciones.extend(filas_bloque)

    return (
        pd.DataFrame(evaluaciones, columns=COLUMNAS_DETECCION),
        pd.DataFrame(candidatos, columns=COLUMNAS_DETECCION),
    )



_detectar_original = detectar_cambios
_cache_local = {}
def detectar_cambios(river_id, theta, ventana_m):
    key = (str(river_id), float(theta), float(ventana_m))
    if key not in _cache_local:
        _cache_local[key] = _detectar_original(*key)
    return tuple(t.copy(deep=True) for t in _cache_local[key])

if EJECUTAR_LOCAL:
    # ── Ejecutar para todos los perfiles cargados ────────────────────────────
    lista_evaluaciones = []
    lista_candidatos = []

    for river_id in perfiles:
        evaluacion, candidatos = detectar_cambios(
            river_id, THETA_DETECCION, VENTANA_LADO_M
        )
        lista_evaluaciones.append(evaluacion)
        lista_candidatos.append(candidatos)

        print(f"{river_id} — {nombres[river_id]}: {len(candidatos)} candidatos")

    evaluaciones_base = pd.concat(lista_evaluaciones, ignore_index=True)
    candidatos_base = pd.concat(lista_candidatos, ignore_index=True)

    # Identificadores válidos dentro de esta ejecución
    candidatos_base["candidato"] = [
        f"R{int(r.ID):02d}_C{j+1:03d}"
        for j, r in candidatos_base.iterrows()
    ]

    display(candidatos_base[
      [
            "candidato", "rio", "tramo", "distancia_m", "elevacion_m",
            "ksn_abajo", "ksn_arriba", "cambio_relativo", "sentido",
        ]
    ])
else:
    evaluaciones_base = pd.DataFrame(columns=COLUMNAS_DETECCION)
    candidatos_base = pd.DataFrame(columns=[*COLUMNAS_DETECCION, 'candidato'])


# Notebook cell 11
if EJECUTAR_SENSIBILIDAD_LOCAL:
    if not EJECUTAR_LOCAL: raise ValueError('Activar primero EJECUTAR_LOCAL.')
    # ── MODIFICAR AQUÍ ──────────────────────────────────────────────────────
    THETAS_COMPARAR = [0.30, 0.45, 0.60]
    VENTANAS_LADO_COMPARAR_M = [250.0, 500.0, 1000.0]

    # Distancia máxima para asociar detecciones de distintas configuraciones
    TOLERANCIA_ASOCIACION_M = 300.0

    # Incluir también la configuración original
    configuraciones = sorted(set(
        [
            (round(float(t), 8), float(w))
            for t in THETAS_COMPARAR
            for w in VENTANAS_LADO_COMPARAR_M
        ]
        + [(THETA_DETECCION, VENTANA_LADO_M)]
    ))

    if any(w <= 0 for _, w in configuraciones):
        raise ValueError("Las ventanas deben ser positivas.")

    for river_id in perfiles:
        for theta, _ in configuraciones:
            if theta not in chi_sensibilidad[river_id]:
                raise ValueError(
                    f"Falta θ={theta} para {river_id}. "
                    "Añádelo al barrido del punto 2."
                )

    listas_pruebas = []
    listas_evaluaciones = []

    for theta, ventana in configuraciones:
        for river_id in perfiles:
            evaluacion, candidatos = detectar_cambios(
                river_id, theta, ventana
            )
            listas_evaluaciones.append(evaluacion)
            listas_pruebas.append(candidatos)

        print(f"Completado: θ={theta:.2f}, ventana por lado={ventana:g} m")

    pruebas_candidatos = pd.concat(listas_pruebas, ignore_index=True)
    evaluaciones_sensibilidad = pd.concat(
        listas_evaluaciones, ignore_index=True
    )

    # ── Asociar cada detección como máximo a un candidato de referencia ──────
    coincidencias = {
        nombre: [] for nombre in candidatos_base["candidato"]
    }
    configuraciones_evaluables = {
        nombre: 0 for nombre in candidatos_base["candidato"]
    }

    for theta, ventana in configuraciones:
        detecciones = pruebas_candidatos[
            pruebas_candidatos["theta"].eq(theta)
            & pruebas_candidatos["ventana_lado_m"].eq(ventana)
        ]

        evaluaciones = evaluaciones_sensibilidad[
            evaluaciones_sensibilidad["theta"].eq(theta)
            & evaluaciones_sensibilidad["ventana_lado_m"].eq(ventana)
        ]

        pares = []

        for _, candidato in candidatos_base.iterrows():
            cerca = evaluaciones[
                evaluaciones["ID"].eq(candidato["ID"])
                & evaluaciones["tramo"].eq(candidato["tramo"])
                & evaluaciones["bloque"].eq(candidato["bloque"])
                & (
                    evaluaciones["distancia_m"] - candidato["distancia_m"]
                ).abs().le(TOLERANCIA_ASOCIACION_M)
            ]

            # Cuenta posiciones evaluadas, aunque el ajuste no superara filtros
            if not cerca.empty:
                configuraciones_evaluables[candidato["candidato"]] += 1

            compatibles = detecciones[
                detecciones["ID"].eq(candidato["ID"])
                & detecciones["tramo"].eq(candidato["tramo"])
                & detecciones["bloque"].eq(candidato["bloque"])
                & detecciones["sentido"].eq(candidato["sentido"])
            ]

            for idx, deteccion in compatibles.iterrows():
                diferencia = abs(
                    deteccion["distancia_m"] - candidato["distancia_m"]
                )
                if diferencia <= TOLERANCIA_ASOCIACION_M:
                    pares.append((
                        diferencia,
                        candidato["candidato"],
                        idx,
                        float(deteccion["distancia_m"]),
                    ))

        usados_base, usados_prueba = set(), set()

        for diferencia, nombre, idx, posicion in sorted(pares):
            if nombre not in usados_base and idx not in usados_prueba:
                coincidencias[nombre].append(posicion)
                usados_base.add(nombre)
                usados_prueba.add(idx)

    catalogo = candidatos_base.copy()

    catalogo["configuraciones_probadas"] = len(configuraciones)
    catalogo["configuraciones_evaluables"] = [
        configuraciones_evaluables[c] for c in catalogo["candidato"]
    ]
    catalogo["detecciones_asociadas"] = [
        len(coincidencias[c]) for c in catalogo["candidato"]
    ]
    catalogo["posicion_min_m"] = [
        min(coincidencias[c]) if coincidencias[c] else np.nan
        for c in catalogo["candidato"]
    ]
    catalogo["posicion_max_m"] = [
        max(coincidencias[c]) if coincidencias[c] else np.nan
        for c in catalogo["candidato"]
    ]

    catalogo["revision"] = "pendiente"
    catalogo["contexto_geologico"] = "no_revisado"
    catalogo["observaciones"] = ""

    display(catalogo[
        [
            "candidato", "rio", "distancia_m",
            "detecciones_asociadas", "configuraciones_probadas",
            "configuraciones_evaluables", "posicion_min_m", "posicion_max_m",
        ]
    ])
else:
    pruebas_candidatos = pd.DataFrame(columns=COLUMNAS_DETECCION)
    evaluaciones_sensibilidad = pd.DataFrame(columns=COLUMNAS_DETECCION)


# Notebook cell 13
if EJECUTAR_TOPOTOOLBOX:
    import numpy as np
    import pandas as pd
    import topotoolbox as tt

    # ── MODIFICAR AQUÍ ────────────────────────────────────────────────
    TOLERANCIA_TT_M = 20.0  # Tolerancia vertical, no cota del knickpoint

    # Mayor tolerancia → normalmente menos candidatos.
    # Valor inicial exploratorio: revisar según la calidad del DEM.


    class Perfil1DTopoToolbox(tt.StreamObject):
        """
        Adapta un perfil ya extraído al algoritmo oficial de TopoToolbox.
        Uso limitado a knickpointfinder, imposemin y lowerenv.
        """

        def __init__(self, distancia):
            d = np.asarray(distancia, dtype=float)

            if (
                len(d) < 3
                or not np.isfinite(d).all()
                or np.any(np.diff(d) <= 0)
            ):
                raise ValueError("Se requieren distancias finitas y crecientes.")

            # Orden interno: cabecera → outlet.
            self.stream = np.arange(len(d), dtype=np.int64)
            self.source = np.arange(len(d) - 1, dtype=np.int64)
            self.target = self.source + 1
            self._distancia = np.ascontiguousarray(
                (d - d[0])[::-1], dtype=np.float32
            )

            if np.any(np.diff(self._distancia) >= 0):
                raise ValueError("Las distancias pierden precisión en float32.")

        def ezgetnal(self, valores, dtype=None):
            a = np.array(valores, dtype=dtype, copy=True)

            if a.shape != self.stream.shape:
                raise ValueError("Elevaciones y nodos no tienen igual longitud.")

            return a

        def upstream_distance(self):
            return self._distancia.copy()

        def node_to_node_distance(self):
            return (
                self._distancia[self.source]
                - self._distancia[self.target]
            )


    if not np.isfinite(TOLERANCIA_TT_M) or TOLERANCIA_TT_M <= 0:
        raise ValueError("TOLERANCIA_TT_M debe ser positiva y finita.")

    filas_tt = []
    filas_qc_tt = []

    # Calcula TODOS los ríos, independientemente de la selección gráfica.
    for river_id, p in perfiles.items():
        x = p["distancia_m"].to_numpy(float)
        z = p["elevacion_m"].to_numpy(float)
        chi = np.asarray(chi_sensibilidad[river_id][THETA_DETECCION])

        for tramo, mascara in mascaras_tramos[river_id].items():
            indices = np.flatnonzero(
                mascara & np.isfinite(z) & np.isfinite(chi)
            )

            # No unir sectores separados por datos faltantes.
            bloques = np.split(
                indices, np.flatnonzero(np.diff(indices) > 1) + 1
            )

            for bloque, idx in enumerate(bloques):
                if len(idx) < 3:
                    continue

                cauce = Perfil1DTopoToolbox(x[idx])
                entrada = np.ascontiguousarray(
                    z[idx][::-1], dtype=np.float32
                )

                z_tt = tt.imposemin(cauce, entrada)[::-1]

                detectados = cauce.knickpointfinder(
                    entrada,
                    knickpoints=np.zeros(len(idx), dtype=bool),
                    tolerance=TOLERANCIA_TT_M,
                )[::-1]

                for k, j in enumerate(idx):
                    filas_qc_tt.append({
                        "ID": river_id,
                        "tramo": tramo,
                        "indice": int(j),
                        "elevacion_original_m": z[j],
                        "elevacion_tt_m": float(z_tt[k]),
                        "correccion_tt_m": float(z_tt[k] - z[j]),
                    })

                    if not detectados[k]:
                        continue

                    filas_tt.append({
                        "ID": river_id,
                        "rio": nombres[river_id],
                        "tramo": tramo,
                        "bloque": bloque,
                        "indice": int(j),
                        "distancia_m": x[j],
                        "elevacion_m": z[j],
                        "chi_m": chi[j],
                        "x_utm": float(p.iloc[j]["x_utm"]),
                        "y_utm": float(p.iloc[j]["y_utm"]),
                        "tolerancia_tt_m": TOLERANCIA_TT_M,
                        "correccion_tt_m": float(z_tt[k] - z[j]),
                    })

    columnas_tt = [
        "ID", "rio", "tramo", "bloque", "indice",
        "distancia_m", "elevacion_m", "chi_m", "x_utm", "y_utm",
        "tolerancia_tt_m", "correccion_tt_m",
    ]

    candidatos_tt = pd.DataFrame(filas_tt, columns=columnas_tt)
    qc_topotoolbox = pd.DataFrame(filas_qc_tt)

    print(f"TopoToolbox: {len(candidatos_tt)} candidatos.")
    print("Las elevaciones originales se conservan en los gráficos.")

    # Descomenta para revisar las correcciones internas:
    display(qc_topotoolbox)
else:
    candidatos_tt = pd.DataFrame(columns=['ID','rio','tramo','bloque','indice','distancia_m',
        'elevacion_m','chi_m','x_utm','y_utm','tolerancia_tt_m','correccion_tt_m'])
    qc_topotoolbox = pd.DataFrame(columns=['ID','tramo','indice','correccion_tt_m'])


# Notebook cell 14
if EJECUTAR_SENSIBILIDAD_TT:
    if not EJECUTAR_TOPOTOOLBOX: raise ValueError('Activar primero EJECUTAR_TOPOTOOLBOX.')
    # ── MODIFICAR AQUÍ ────────────────────────────────────────────────
    TOLERANCIAS_TT_M = [10.0, 20.0, 40.0]

    # Distancia sobre el cauce para asociar candidatos entre ejecuciones.
    # Es un parámetro exploratorio: revisar según resolución y morfología.
    DISTANCIA_ASOCIACION_TT_M = 300.0


    def ejecutar_tt_tolerancia(tolerancia):
        """Ejecuta el método oficial sobre los mismos perfiles y tramos."""
        filas = []

        for river_id, p in perfiles.items():
            x = p["distancia_m"].to_numpy(float)
            z = p["elevacion_m"].to_numpy(float)
            chi = np.asarray(
                chi_sensibilidad[river_id][THETA_DETECCION]
            )

            for tramo, mascara in mascaras_tramos[river_id].items():
                indices = np.flatnonzero(
                    mascara & np.isfinite(z) & np.isfinite(chi)
                )
                bloques = np.split(
                    indices,
                    np.flatnonzero(np.diff(indices) > 1) + 1,
                )

                for bloque, idx in enumerate(bloques):
                    if len(idx) < 3:
                        continue

                    cauce = Perfil1DTopoToolbox(x[idx])
                    entrada = np.ascontiguousarray(
                        z[idx][::-1], dtype=np.float32
                    )

                    detectados = cauce.knickpointfinder(
                        entrada,
                        knickpoints=np.zeros(len(idx), dtype=bool),
                        tolerance=float(tolerancia),
                    )[::-1]

                    for j in idx[detectados]:
                        filas.append({
                            "ID": str(river_id),
                            "tramo": tramo,
                            "bloque": bloque,
                            "indice": int(j),
                            "distancia_m": float(x[j]),
                            "tolerancia_tt_m": float(tolerancia),
                        })

        return pd.DataFrame(filas, columns=[
            "ID", "tramo", "bloque", "indice",
            "distancia_m", "tolerancia_tt_m",
        ])


    tolerancias = sorted(set(
        map(float, TOLERANCIAS_TT_M + [TOLERANCIA_TT_M])
    ))

    if any(not np.isfinite(t) or t <= 0 for t in tolerancias):
        raise ValueError("Las tolerancias deben ser positivas y finitas.")

    if (
        not np.isfinite(DISTANCIA_ASOCIACION_TT_M)
        or DISTANCIA_ASOCIACION_TT_M < 0
    ):
        raise ValueError("Revisa DISTANCIA_ASOCIACION_TT_M.")

    resultados_tt = []

    for tolerancia in tolerancias:
        # Reutilizar la ejecución de referencia.
        if np.isclose(tolerancia, TOLERANCIA_TT_M):
            tabla = candidatos_tt[
                ["ID", "tramo", "bloque", "indice", "distancia_m"]
            ].copy()
            tabla["ID"] = tabla["ID"].astype(str)
            tabla["tolerancia_tt_m"] = tolerancia
        else:
            tabla = ejecutar_tt_tolerancia(tolerancia)

        resultados_tt.append(tabla)
        print(f"Tolerancia {tolerancia:g} m: {len(tabla)} candidatos.")

    pruebas_tt = pd.concat(resultados_tt, ignore_index=True)

    # Asociaciones uno a uno por configuración, priorizando menor distancia.
    referencia_tt = candidatos_tt.reset_index(drop=True).copy()
    referencia_tt["ID"] = referencia_tt["ID"].astype(str)

    asociaciones = {i: [] for i in referencia_tt.index}

    for tolerancia in tolerancias:
        prueba = pruebas_tt[
            pruebas_tt["tolerancia_tt_m"].eq(tolerancia)
        ]

        pares = []

        for i, base in referencia_tt.iterrows():
            compatibles = prueba[
                prueba["ID"].eq(base["ID"])
                & prueba["tramo"].eq(base["tramo"])
                & prueba["bloque"].eq(base["bloque"])
            ]

            for j, candidato in compatibles.iterrows():
                diferencia = abs(
                    candidato["distancia_m"] - base["distancia_m"]
                )

                if diferencia <= DISTANCIA_ASOCIACION_TT_M:
                    pares.append((
                        diferencia, i, j, candidato["distancia_m"]
                    ))

        usados_base, usados_prueba = set(), set()

        for diferencia, i, j, posicion in sorted(pares):
            if i in usados_base or j in usados_prueba:
                continue

            asociaciones[i].append((tolerancia, posicion))
            usados_base.add(i)
            usados_prueba.add(j)

    sensibilidad_tt = referencia_tt.copy()

    sensibilidad_tt["n_tolerancias"] = len(tolerancias)
    sensibilidad_tt["n_detecciones"] = [
        len(asociaciones[i]) for i in referencia_tt.index
    ]
    sensibilidad_tt["tolerancias_detectado"] = [
        ", ".join(f"{t:g}" for t, _ in asociaciones[i])
        for i in referencia_tt.index
    ]
    sensibilidad_tt["posicion_min_m"] = [
        min(p for _, p in asociaciones[i])
        if asociaciones[i] else np.nan
        for i in referencia_tt.index
    ]
    sensibilidad_tt["posicion_max_m"] = [
        max(p for _, p in asociaciones[i])
        if asociaciones[i] else np.nan
        for i in referencia_tt.index
    ]

    display(sensibilidad_tt[[
        "ID", "tramo", "indice", "distancia_m",
        "n_detecciones", "n_tolerancias",
        "tolerancias_detectado", "posicion_min_m", "posicion_max_m",
    ]])

    # La frecuencia de detección NO es una probabilidad de que sea real.
else:
    sensibilidad_tt = pd.DataFrame(columns=['ID','tramo','indice','n_detecciones','n_tolerancias',
        'tolerancias_detectado','posicion_min_m','posicion_max_m'])


# Notebook cell 16
from pathlib import Path
from scipy.spatial import cKDTree

# ── MODIFICAR CUANDO EXISTAN RESULTADOS DE LSDTopoTools ─────────────
# Los archivos deben estar en el mismo CRS que tus perfiles: EPSG:32718.
# Usa la misma concavidad que THETA_DETECCION.
#
# Ejemplo de formato; los identificadores son solo ilustrativos:
#
# RESULTADOS_GAILLETON = {
#     "6": {
#         "archivo": r"C:\ruta\Reloca_ksnkp.csv",
#         "source_key": 0,
#         "basin_key": 0,
#         "theta": 0.45,
#         "crs": "EPSG:32718",
#     }
# }

RESULTADOS_GAILLETON = {}

# Distancia máxima entre un candidato de LSDTopoTools y tu perfil.
# Revisar según resolución del DEM y diferencias entre redes.
DISTANCIA_MAX_G_M = 30.0

filas_g = []
rechazos_g = []
estado_gailleton = {}

for river_id, p in perfiles.items():
    if river_id not in RESULTADOS_GAILLETON:
        estado_gailleton[river_id] = "pendiente"
        continue

    cfg = RESULTADOS_GAILLETON[river_id]

    if cfg["crs"] != "EPSG:32718":
        raise ValueError("Transforma y verifica las coordenadas antes de importar.")

    if not np.isclose(cfg["theta"], THETA_DETECCION):
        raise ValueError(f"Río {river_id}: las concavidades no coinciden.")

    archivo = Path(cfg["archivo"])

    if not archivo.name.endswith("_ksnkp.csv"):
        raise ValueError("Selecciona el archivo final terminado en _ksnkp.csv.")

    datos = pd.read_csv(archivo)
    requeridas = {
        "X", "Y", "elevation", "source_key", "basin_key",
        "delta_ksn", "delta_segelev",
    }

    if requeridas - set(datos.columns):
        raise ValueError(f"Faltan columnas: {requeridas - set(datos.columns)}")

    # Verificar que el cauce fue procesado, incluso si no tuvo candidatos.
    archivo_segmentos = archivo.with_name(
        archivo.name.replace("_ksnkp.csv", "_ksnkp_mchi.csv")
    )
    segmentos = pd.read_csv(archivo_segmentos)

    canal = (
        segmentos["source_key"].eq(cfg["source_key"])
        & segmentos["basin_key"].eq(cfg["basin_key"])
    )

    if not canal.any():
        raise ValueError(f"Río {river_id}: ese cauce no aparece en los resultados.")

    datos = datos[
        datos["source_key"].eq(cfg["source_key"])
        & datos["basin_key"].eq(cfg["basin_key"])
    ]

    xy = p[["x_utm", "y_utm"]].to_numpy(float)

    if not np.isfinite(xy).all():
        raise ValueError(f"Río {river_id}: coordenadas inválidas.")

    arbol = cKDTree(xy)

    for fila, r in datos.iterrows():
        registro = {
            "ID": river_id,
            "archivo": str(archivo),
            "fila_lsd": int(fila),
            "x_lsd": r["X"],
            "y_lsd": r["Y"],
            "elevacion_lsd_m": r["elevation"],
            "delta_ksn_lsd": r["delta_ksn"],
            "delta_segelev_lsd": r["delta_segelev"],
        }

        if not np.isfinite([r["X"], r["Y"]]).all():
            rechazos_g.append({
                **registro, "motivo": "coordenadas_invalidas"
            })
            continue

        distancia, j = arbol.query([r["X"], r["Y"]])
        j = int(j)

        tramos = [
            t for t, mascara in mascaras_tramos[river_id].items()
            if mascara[j]
        ]

        if (
            distancia > DISTANCIA_MAX_G_M
            or len(tramos) != 1
            or not np.isfinite(p.iloc[j]["elevacion_m"])
        ):
            rechazos_g.append({
                **registro,
                "desajuste_xy_m": distancia,
                "motivo": "sin_asociacion_valida_al_perfil",
            })
            continue

        punto = p.iloc[j]

        filas_g.append({
            **registro,
            "rio": nombres[river_id],
            "tramo": tramos[0],
            "indice": j,
            "distancia_m": punto["distancia_m"],
            "elevacion_m": punto["elevacion_m"],
            "chi_m": chi_sensibilidad[river_id][THETA_DETECCION][j],
            "x_utm": punto["x_utm"],
            "y_utm": punto["y_utm"],
            "desajuste_xy_m": distancia,
        })

    estado_gailleton[river_id] = "importado"

columnas_g = [
    "ID", "rio", "tramo", "indice", "distancia_m", "elevacion_m",
    "chi_m", "x_utm", "y_utm", "archivo", "fila_lsd",
    "x_lsd", "y_lsd", "elevacion_lsd_m",
    "delta_ksn_lsd", "delta_segelev_lsd", "desajuste_xy_m",
]

candidatos_gailleton = pd.DataFrame(filas_g, columns=columnas_g)
rechazos_gailleton = pd.DataFrame(rechazos_g)

print(f"Gailleton: {len(candidatos_gailleton)} candidatos importados.")
print(
    "Ríos pendientes:",
    sum(estado == "pendiente" for estado in estado_gailleton.values()),
)


# Notebook cell 18
tablas = []
for metodo, tabla, prefijo in [('Local', candidatos_base, 'L'), ('TopoToolbox', candidatos_tt, 'T'),
                               ('Gailleton', candidatos_gailleton, 'G')]:
    t = tabla.copy(); t['metodo'] = metodo
    t['perfil_hash'] = t['ID'].astype(str).map(lambda i: pa.profile_hash(perfiles[i]))
    t['id_kp'] = [f"R{int(r.ID):02d}_{prefijo}_I{int(r.indice):05d}_{r.perfil_hash[:8]}"
                  + (f'_{j:04d}' if metodo == 'Gailleton' else '') for j,r in enumerate(t.itertuples())]
    tablas.append(t)
manuales, manuales_rechazados = pa.load_manuals(perfiles)
manuales['rio'] = manuales.ID.map(nombres)
manuales['chi_m'] = [chi_sensibilidad[r.ID][THETA_DETECCION][int(r.indice)] for r in manuales.itertuples()]
manuales['tramo'] = [next((t for t,m in mascaras_tramos[r.ID].items() if m[int(r.indice)]), 'fuera_de_tramos') for r in manuales.itertuples()]
catalogo_revision = pd.concat([*tablas, manuales], ignore_index=True)
catalogo_revision['ID'] = catalogo_revision.ID.astype(str)
catalogo_revision['me_interesa'] = catalogo_revision.metodo.eq('Manual')
catalogo_revision['revision'] = 'pendiente'
catalogo_revision['observaciones'] = ''
if catalogo_revision.id_kp.duplicated().any(): raise ValueError('IDs de KP duplicados.')

# Bloques espaciales comunes: no asociar a través de huecos de elevación.
bloques = {}
for rid,p in perfiles.items():
    fin = np.isfinite(p.elevacion_m.to_numpy(float))
    bloques[rid] = np.cumsum(fin & ~np.r_[False, fin[:-1]])
catalogo_revision['bloque_comun'] = [int(bloques[r.ID][int(r.indice)]) for r in catalogo_revision.itertuples()]

if CFG['revision_candidatos_csv']:
    ruta_revision = Path(CFG['revision_candidatos_csv'])
    decisiones = pd.read_csv(ruta_revision, dtype={'ID': str})
    if decisiones.id_kp.duplicated().any(): raise ValueError('Revisión con IDs duplicados.')
    for row in decisiones.to_dict('records'):
        match = catalogo_revision.id_kp.eq(row['id_kp']) & catalogo_revision.perfil_hash.eq(row['perfil_hash'])
        for col in ['me_interesa','revision','observaciones']:
            if col in row: catalogo_revision.loc[match,col] = row[col]
    deps[str(ruta_revision)] = pa.file_hash(ruta_revision)
display(manuales_rechazados)


# Notebook cell 20
from itertools import combinations

# ── MODIFICAR AQUÍ ────────────────────────────────────────────────
DISTANCIA_COINCIDENCIA_M = 300.0

if (
    not np.isfinite(DISTANCIA_COINCIDENCIA_M)
    or DISTANCIA_COINCIDENCIA_M < 0
):
    raise ValueError("La distancia debe ser finita y no negativa.")

filas_coincidencias = []

for (river_id, tramo, bloque), grupo in catalogo_revision.groupby(
    ["ID", "tramo", "bloque_comun"], sort=False
):
    metodos = sorted(grupo["metodo"].unique())

    for metodo_a, metodo_b in combinations(metodos, 2):
        a = grupo[grupo["metodo"].eq(metodo_a)]
        b = grupo[grupo["metodo"].eq(metodo_b)]

        posibles = []

        for _, ca in a.iterrows():
            for _, cb in b.iterrows():
                separacion = abs(
                    ca["distancia_m"] - cb["distancia_m"]
                )

                if separacion <= DISTANCIA_COINCIDENCIA_M:
                    posibles.append((
                        separacion, ca["id_kp"], cb["id_kp"]
                    ))

        # Un candidato se asocia como máximo a uno del otro método.
        usados_a, usados_b = set(), set()

        for separacion, id_a, id_b in sorted(posibles):
            if id_a in usados_a or id_b in usados_b:
                continue

            filas_coincidencias.append({
                "ID": river_id,
                "tramo": tramo,
                "metodo_a": metodo_a,
                "id_a": id_a,
                "metodo_b": metodo_b,
                "id_b": id_b,
                "separacion_m": separacion,
            })

            usados_a.add(id_a)
            usados_b.add(id_b)

coincidencias_metodos = pd.DataFrame(
    filas_coincidencias,
    columns=[
        "ID", "tramo", "metodo_a", "id_a",
        "metodo_b", "id_b", "separacion_m",
    ],
)

# Guardar las asociaciones directas sin fusionar candidatos.
vecinos = {
    id_kp: [] for id_kp in catalogo_revision["id_kp"]
}

for _, fila in coincidencias_metodos.iterrows():
    vecinos[fila["id_a"]].append(fila["id_b"])
    vecinos[fila["id_b"]].append(fila["id_a"])

catalogo_revision["coincide_con"] = [
    "; ".join(vecinos[id_kp])
    for id_kp in catalogo_revision["id_kp"]
]

display(coincidencias_metodos)

# Estas asociaciones no agrupan automáticamente tres puntos
# como un único knickpoint ni indican un origen geológico.


# Notebook cell 22
revision_no_aplicada = pd.DataFrame(columns=['id_kp','perfil_hash','motivo'])
if CFG['revision_candidatos_csv']:
    pares_vigentes = set(zip(catalogo_revision.id_kp, catalogo_revision.perfil_hash))
    revision_no_aplicada = decisiones.loc[[
        (r.id_kp,r.perfil_hash) not in pares_vigentes for r in decisiones.itertuples()]].copy()
    revision_no_aplicada['motivo'] = 'ID o huella sin correspondencia en esta corrida'


# Notebook cell 24
estados = pd.DataFrame([{'ID': rid, 'rio': nombres[rid],
    'manuales_validos': int(manuales.ID.eq(rid).sum()),
    'manuales_incompatibles': int((manuales_rechazados.ID.astype(str).eq(rid) & manuales_rechazados.estado.eq('incompatible')).sum()),
    'estado_manual': 'ok' if manuales.ID.eq(rid).any() else 'sin_manual_valido',
    'local': 'ejecutado' if EJECUTAR_LOCAL else 'no_ejecutado',
    'topotoolbox': 'ejecutado' if EJECUTAR_TOPOTOOLBOX else 'no_ejecutado',
    'gailleton': estado_gailleton.get(rid, 'no_ejecutado'),
    'sensibilidad_local': 'ejecutada' if EJECUTAR_SENSIBILIDAD_LOCAL else 'no_ejecutada',
    'sensibilidad_tt': 'ejecutada' if EJECUTAR_SENSIBILIDAD_TT else 'no_ejecutada'} for rid in pa.EXPECTED])
tablas_exportar = {'catalogo_candidatos': catalogo_revision, 'manuales': manuales,
    'manuales_rechazados': manuales_rechazados, 'estados': estados, 'qc': control_calidad,
    'ajustes': ajustes_segmentos, 'barrido_theta': tabla_sensibilidad,
    'coincidencias': coincidencias_metodos, 'sensibilidad_tt': sensibilidad_tt,
    'sensibilidad_local': pruebas_candidatos, 'evaluaciones_sensibilidad_local': evaluaciones_sensibilidad,
    'evaluaciones_local': evaluaciones_base, 'qc_topotoolbox': qc_topotoolbox}
tablas_exportar['revision_no_aplicada'] = revision_no_aplicada
archivos = {k: pa.save_table('analisis_perfiles', k, t) for k,t in tablas_exportar.items()}
pa.finish('analisis_perfiles', archivos, deps, extra={'tramos': TRAMOS,
    'theta_ajuste': THETA_AJUSTE, 'A0_m2': A0_ANALISIS,
    'detector_local': {'ventana_lado_m': VENTANA_LADO_M, 'paso_busqueda_m': PASO_BUSQUEDA_M,
        'separacion_m': SEPARACION_CANDIDATOS_M, 'min_puntos': MIN_PUNTOS_LADO,
        'r2_min': R2_MINIMO, 'cambio_min': CAMBIO_MINIMO, 'prominencia_min': PROMINENCIA_MINIMA}})
display(estados)


# Notebook cell 26
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ── MODIFICAR AQUÍ: SOLO VISUALIZACIÓN ─────────────────────────────
#RIOS_VER = ["6"]
# RIOS_VER = ["1", "3", "16"]
RIOS_VER = ['6']

METODOS_VER = ["Local", "TopoToolbox", "Gailleton"]

MOSTRAR_ETIQUETAS = True
MOSTRAR_FRANJA = True
ELEVACION_MIN = 180
ELEVACION_MAX = 220

ESTILOS = {
    "Local":      {"color": "#D55E00", "marker": "o", "s": 65},
    "TopoToolbox": {"color": "#7B2CBF", "marker": "^", "s": 100},
    "Gailleton":  {"color": "#009E73", "marker": "s", "s": 135},
}

ESTILOS['Manual'] = {'color': 'black', 'marker': '*', 's': 150}
METODOS_VER.append('Manual')
catalogo_metodos = catalogo_revision.copy()
catalogo_metodos['etiqueta'] = catalogo_metodos['id_kp']
def figura_revision(river_id):
    river_id = str(river_id)

    if river_id not in perfiles:
        raise ValueError(f"No está cargado el río {river_id}.")

    if set(METODOS_VER) - set(ESTILOS):
        raise ValueError("Revisa los nombres en METODOS_VER.")

    p = perfiles[river_id]
    chi = chi_sensibilidad[river_id][THETA_DETECCION]
    candidatos = catalogo_metodos[
        catalogo_metodos["ID"].eq(river_id)
    ]

    fig, axs = plt.subplots(1, 3, figsize=(18, 5.7))

    axs[0].plot(
        p["distancia_m"] / 1000, p["elevacion_m"],
        color="steelblue", linewidth=1.2,
    )
    axs[1].plot(
        chi, p["elevacion_m"],
        color="steelblue", linewidth=1.2,
    )
    axs[2].plot(
        p["x_utm"], p["y_utm"],
        color="steelblue", linewidth=1.2,
    )
    axs[2].scatter(
        p.iloc[0]["x_utm"], p.iloc[0]["y_utm"],
        marker="*", s=100, color="black", zorder=6,
    )

    leyenda = []

    for orden, metodo in enumerate(METODOS_VER):
        estilo = ESTILOS[metodo]
        grupo = candidatos[candidatos["metodo"].eq(metodo)]

        pendiente = (
            metodo == "Gailleton"
            and estado_gailleton[river_id] == "pendiente"
        )
        cantidad = "pendiente" if pendiente else str(len(grupo))

        leyenda.append(Line2D(
            [], [],
            linestyle="none",
            marker=estilo["marker"],
            color=estilo["color"],
            markerfacecolor="none",
            markersize=8,
            label=f"{metodo}: {cantidad}",
        ))

        for _, candidato in grupo.iterrows():
            puntos = [
                (
                    candidato["distancia_m"] / 1000,
                    candidato["elevacion_m"],
                ),
                (candidato["chi_m"], candidato["elevacion_m"]),
                (candidato["x_utm"], candidato["y_utm"]),
            ]

            for ax, (x, y) in zip(axs, puntos):
                # Símbolos vacíos para distinguir métodos coincidentes.
                ax.scatter(
                    x, y,
                    marker=estilo["marker"],
                    s=estilo["s"],
                    facecolors="none",
                    edgecolors=estilo["color"],
                    linewidths=1.6,
                    zorder=7 + orden,
                )

                if MOSTRAR_ETIQUETAS:
                    ax.annotate(
                        candidato["etiqueta"],
                        (x, y),
                        xytext=(5, 8 + 10 * orden),
                        textcoords="offset points",
                        fontsize=7,
                        color=estilo["color"],
                    )

    if MOSTRAR_FRANJA:
        # Solo perfiles: el mapa tiene Norte, no elevación, en el eje Y.
        for ax in axs[:2]:
            ax.axhspan(
                ELEVACION_MIN, ELEVACION_MAX,
                color="gold", alpha=0.22, zorder=0,
            )
            for elevacion in (ELEVACION_MIN, ELEVACION_MAX):
                ax.axhline(
                    elevacion, color="darkorange",
                    linestyle="--", linewidth=0.8,
                )

    axs[0].set(
        title="Perfil longitudinal",
        xlabel="Distancia desde el outlet (km)",
        ylabel="Elevación (m)",
    )
    axs[1].set(
        title=f"χ–elevación | θ={THETA_DETECCION}",
        xlabel="χ (m)",
        ylabel="Elevación (m)",
    )
    axs[2].set(
        title="Ubicación sobre el cauce",
        xlabel="Este (m)",
        ylabel="Norte (m)",
    )
    axs[2].set_aspect("equal")
    axs[2].ticklabel_format(style="plain", useOffset=False)

    for ax in axs:
        ax.grid(alpha=0.25)

    fig.suptitle(
        f"{river_id} — {nombres[river_id]} | "
        "Candidatos por método"
    )
    fig.legend(
        handles=leyenda,
        loc="lower center",
        ncol=len(METODOS_VER),
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))

    return fig


for river_id in map(str, RIOS_VER):
    fig = figura_revision(river_id)
    plt.show()
    plt.close(fig)


# Notebook cell 28
EDITAR_MANUALES = False
RIOS_EDITAR = ['6']
if EDITAR_MANUALES:
    from datetime import datetime, timezone
    get_ipython().run_line_magic('matplotlib','tk')
    propuestas = manuales.to_dict('records')
    try:
        for rid in map(str, RIOS_EDITAR):
            p = perfiles[rid]; chi = chi_sensibilidad[rid][THETA_DETECCION]
            z = p.elevacion_m.to_numpy(float)
            fig,ax = plt.subplots(figsize=(12,5))
            ax.plot(chi,z);ax.set(title=nombres[rid]+' — clics para añadir; Enter termina',xlabel='χ (m)',ylabel='Elevación (m)')
            fig.canvas.draw()
            clics = plt.ginput(n=-1,timeout=0)
            validos = np.flatnonzero(np.isfinite(z) & np.isfinite(chi))
            puntos_pantalla = ax.transData.transform(np.column_stack([chi[validos],z[validos]]))
            for xc,zc in clics:
                pixel = ax.transData.transform((xc,zc))
                j = int(validos[np.argmin(np.linalg.norm(puntos_pantalla-pixel,axis=1))])
                h = pa.profile_hash(p); id_kp = f'R{int(rid):02d}_M_I{j:05d}_{h[:8]}'
                if any(r['id_kp']==id_kp for r in propuestas): continue
                row=p.iloc[j]
                propuestas.append({'ID':rid,'id_kp':id_kp,'indice':j,'perfil_hash':h,
                    'distancia_m':float(row.distancia_m),'elevacion_m':float(row.elevacion_m),
                    'x_utm':float(row.x_utm),'y_utm':float(row.y_utm),
                    'fecha_seleccion':datetime.now(timezone.utc).isoformat(),
                    'theta_seleccion':THETA_DETECCION,'A0_seleccion_m2':A0_ANALISIS})
            plt.close(fig)
        fecha_edicion = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        ruta_propuesta = pa.save_table('analisis_perfiles','manuales_propuestos_'+fecha_edicion,pd.DataFrame(propuestas))
        print('Propuesta guardada:', ruta_propuesta)
        print('Para adoptarla: kp_manual_csv = esa ruta; kp_separador = coma; kp_decimal = punto.')
    finally:
        plt.close('all')
        get_ipython().run_line_magic('matplotlib','inline')


# Notebook cell 29
# CELDA PARA PEGAR EN JUPYTER. Preparada sin ejecutar cálculos.
# Ejecutar después de cargar perfiles; no requiere KP ni volver a cargar el DEM.
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.stats import linregress
from IPython.display import display

# ===================== MODIFICAR AQUÍ =====================
THETA_K = 0.40
A0_K_M2 = 1e6
MIN_PUNTOS_K = 10
GUARDAR_TABLA = False

# Completar cuando se confirme el ranking espacial. No se fijan cinco por suposición.
RIOS_K = ["3", "4", "7", "8", "9"]  # Ejemplo de formato: ["ID", "ID", "ID", "ID", "ID"]

# Asociaciones recuperadas de las tablas: CANDIDATAS, no cinco ríos aprobados.
# U_m_ka se deja vacía: la tasa del marcador fuente no se transfiere automáticamente.
CONFIG_K = {
    "3": {"terraza_id": "JM2015_028", "cota_terraza_m": 67.80, "muestra": "CONT 11", "OD_pct": 11.70},
    "4": {"terraza_id": "JM2015_030", "cota_terraza_m": 85.60, "muestra": "CONT 11", "OD_pct": 11.70},
    "7": {"terraza_id": "JM2015_045", "cota_terraza_m": 62.68, "muestra": "PUHE 11", "OD_pct": 8.98},
    "8": {"terraza_id": "JM2015_046", "cota_terraza_m": 93.47, "muestra": "PUHE 11", "OD_pct": 8.98},
    "9": {"terraza_id": "JM2015_048", "cota_terraza_m": 70.74, "muestra": "PUHE 11", "OD_pct": 8.98},
}
for _cfg in CONFIG_K.values():
    _cfg.update(
        asociacion_confirmada=False,
        edad_adoptada_ka=None,  # Distinguir edad medida de edad asignada/correlacionada.
        fuente_edad="",
        U_m_ka=None,
        fuente_U="",  # Documentar también tratamiento del nivel marino, si corresponde.
        equilibrio_adoptado=False,
        cruce_id=None,  # Si hay varios cruces por cota, elegir uno de la tabla que se mostrará.
        fin_planicie_m=None,  # Distancia desde TU outlet hasta el inicio de roca confirmado.
        evidencia_planicie="",  # Capa, observación y criterio; pendiente baja sola no basta.
    )

# Escribe tus decisiones DESPUÉS del bucle anterior. Ejemplo de formato:
CONFIG_K["3"].update(asociacion_confirmada=True, edad_adoptada_ka=..., fuente_edad="...",
     U_m_ka=..., fuente_U="...", equilibrio_adoptado=True,
     cruce_id=..., fin_planicie_m=..., evidencia_planicie="...")
fin_planicie_m = 0
#fin_planicie_m = None: evaluación pendiente; NO implica ausencia de planicie.
# Esta celda recorta el tramo inferior contiguo al outlet. No une tramos separados
# por planicies interiores; esos casos requieren definir ajustes independientes.
# =========================================================

def _cruces_cota_K(p, cota):
    s = p["distancia_m"].to_numpy(float)
    z = p["elevacion_m"].to_numpy(float)
    puntos, planos = [], []
    for j in range(len(s) - 1):
        if not np.isfinite(z[j:j + 2]).all():
            continue
        if z[j] == z[j + 1]:
            if z[j] == cota:
                planos.append((s[j], s[j + 1]))
        elif min(z[j], z[j + 1]) <= cota <= max(z[j], z[j + 1]):
            f = (cota - z[j]) / (z[j + 1] - z[j])
            puntos.append(s[j] + f * (s[j + 1] - s[j]))
    return np.unique(np.round(puntos, 6)), planos


def _ajustar_K(p, inicio, fin, U_m_ka):
    s, z, area = (p[c].to_numpy(float) for c in ["distancia_m", "elevacion_m", "area_m2"])
    if len(p) < 2 or not np.isfinite(s).all() or np.any(np.diff(s) <= 0) or not np.isclose(s[0], 0):
        raise ValueError("Distancia inválida o el perfil no empieza en tu outlet.")
    if not np.isfinite(area).all() or np.any(area <= 0):
        raise ValueError("Área inválida; no eliminar filas para cerrar huecos.")
    if not np.isfinite([inicio, fin]).all() or not s[0] <= inicio < fin <= s[-1]:
        raise ValueError("Límites fuera del perfil o intervalo vacío.")
    lo = max(0, np.searchsorted(s, inicio, side="right") - 1)
    hi = min(len(s) - 1, np.searchsorted(s, fin))
    if not np.isfinite(z[lo:hi + 1]).all():
        raise ValueError("El intervalo contiene huecos de elevación.")
    nodos = (s >= inicio) & (s <= fin)
    if int(nodos.sum()) < MIN_PUNTOS_K:
        raise ValueError(f"Menos de {MIN_PUNTOS_K} nodos originales en el intervalo.")

    # Recalcular solo χ 1D con theta=0.40; no reutilizar chi_m si era theta=0.45.
    chi = cumulative_trapezoid((A0_K_M2 / area) ** THETA_K, x=s, initial=0)
    posiciones = np.unique(np.r_[inicio, s[nodos], fin])
    x = np.interp(posiciones, s, chi)
    x -= x[0]
    alturas = np.interp(posiciones, s, z)
    fit = linregress(x, alturas)
    if not np.isfinite(fit.slope) or fit.slope <= 0:
        raise ValueError("Pendiente z–χ no positiva; no se utiliza su valor absoluto.")
    ksn = fit.slope * A0_K_M2 ** THETA_K
    return {
        "inicio_m": inicio, "fin_m": fin, "longitud_m": fin - inicio,
        "n_nodos_originales": int(nodos.sum()), "M_chi": fit.slope,
        "ksn_calibracion": ksn, "R2": fit.rvalue ** 2,
        "K": (U_m_ka / 1000) / ksn,
        "unidad_K": f"m^{1 - 2 * THETA_K:g}/año",
    }


if not RIOS_K:
    print("Falta confirmar los cinco IDs por proximidad. Estas asociaciones son candidatas:")
    display(pd.DataFrame.from_dict(CONFIG_K, orient="index").rename_axis("ID").reset_index())
else:
    ids_K = list(map(str, RIOS_K))
    if len(ids_K) != 5 or len(set(ids_K)) != 5:
        raise ValueError("RIOS_K debe contener cinco IDs diferentes confirmados.")
    if "perfiles_ksn" in globals():
        perfiles_para_K = perfiles_ksn
    elif "perfiles" in globals():
        perfiles_para_K = perfiles
    else:
        raise NameError("Cargar primero los perfiles exportados; no hace falta repetir hidrología.")

    filas_K, filas_cruces_K = [], []
    for rid in ids_K:
        cfg = CONFIG_K.get(rid, {})
        base = {"ID": rid, "terraza_id": cfg.get("terraza_id", ""),
                "muestra": cfg.get("muestra", ""), "OD_pct": cfg.get("OD_pct"),
                "edad_adoptada_ka": cfg.get("edad_adoptada_ka"),
                "fuente_edad": cfg.get("fuente_edad", ""),
                "U_m_ka": cfg.get("U_m_ka"), "fuente_U": cfg.get("fuente_U", ""),
                "theta": THETA_K, "A0_m2": A0_K_M2, "n": 1, "K": np.nan}
        try:
            p = perfiles_para_K[rid]
            puntos, planos = _cruces_cota_K(p, float(cfg["cota_terraza_m"]))
            for i, posicion in enumerate(puntos, 1):
                filas_cruces_K.append({"ID": rid, "cruce_id": i, "distancia_m": posicion})
            if not cfg.get("asociacion_confirmada"):
                raise ValueError("Asociación muestra–terraza–río pendiente de confirmar.")
            od = float(cfg["OD_pct"])
            if not np.isfinite(od) or not 0 <= od < 15:
                raise ValueError("La muestra no cumple OD estrictamente menor a 15%.")
            edad, U = float(cfg["edad_adoptada_ka"]), float(cfg["U_m_ka"])
            if not np.isfinite([edad, U]).all() or edad <= 0 or U <= 0:
                raise ValueError("Faltan edad adoptada y tasa U positivas.")
            if not cfg.get("fuente_edad", "").strip() or not cfg.get("fuente_U", "").strip():
                raise ValueError("Falta documentar la edad y la tasa; no se calcula U desde la edad sola.")
            if not cfg.get("equilibrio_adoptado"):
                raise ValueError("Falta adoptar el supuesto de equilibrio para calibrar K.")
            if planos:
                raise ValueError("La cota coincide con un segmento plano; revisar límite.")
            if len(puntos) == 0:
                raise ValueError("No hay cruce con la cota de terraza.")
            elegido = cfg.get("cruce_id")
            if elegido is None:
                if len(puntos) != 1:
                    raise ValueError("Hay varios cruces; elegir cruce_id en la tabla.")
                elegido = 1
            elegido = float(elegido)
            if not elegido.is_integer() or not 1 <= elegido <= len(puntos):
                raise ValueError("cruce_id inválido.")
            fin = float(puntos[int(elegido) - 1])
            for escenario, inicio in [("sin_recorte", 0.0), ("con_recorte_aluvial", cfg.get("fin_planicie_m"))]:
                row = {**base, "escenario": escenario, "cruce_id": int(elegido),
                       "estado": "pendiente", "detalle": "", "evidencia_planicie": cfg.get("evidencia_planicie", "")}
                try:
                    if inicio is None:
                        raise ValueError("Planicie no evaluada: escenario recortado pendiente.")
                    if escenario == "con_recorte_aluvial" and not cfg.get("evidencia_planicie", "").strip():
                        raise ValueError("Falta evidencia de la exclusión o de ausencia de tramo aluvial inferior.")
                    row.update(_ajustar_K(p, float(inicio), fin, U))
                    row.update(estado="K_condicional", detalle="n=1; equilibrio adoptado; verificar aplicabilidad geológica")
                    if escenario == "con_recorte_aluvial" and float(inicio) == 0:
                        row["detalle"] += "; sin tramo inferior que excluir: mismo intervalo"
                except (ValueError, TypeError) as exc:
                    row["detalle"] = str(exc)
                filas_K.append(row)
        except (KeyError, ValueError, TypeError) as exc:
            for escenario in ["sin_recorte", "con_recorte_aluvial"]:
                filas_K.append({**base, "escenario": escenario, "estado": "pendiente", "detalle": str(exc)})

    cruces_K = pd.DataFrame(filas_cruces_K, columns=["ID", "cruce_id", "distancia_m"])
    comparacion_K = pd.DataFrame(filas_K)
    display(cruces_K)
    display(comparacion_K)
    if GUARDAR_TABLA:
        fecha_K = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        carpeta_K = Path.cwd() / "resultados_K" / fecha_K
        carpeta_K.mkdir(parents=True, exist_ok=True)
        comparacion_K.to_csv(carpeta_K / f"comparacion_K_{fecha_K}.csv", index=False, encoding="utf-8-sig")
        cruces_K.to_csv(carpeta_K / f"cruces_K_{fecha_K}.csv", index=False, encoding="utf-8-sig")
        print("Guardado en:", carpeta_K)


# Notebook cell 30


# Notebook cell 31

