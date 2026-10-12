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

pa.NOTEBOOK_PATH = CARPETA_NOTEBOOKS / 'AnalisisTerrazas_20261005.ipynb'


# Notebook cell 3
perfiles, nombres, manifest = pa.load_profiles()
manuales, manuales_rechazados = pa.load_manuals(perfiles)
deps = pa.dependency('perfiles')
for key in ['kp_manual_csv', 'terrazas_excel', 'revision_muestras_excel', 'seleccion_excel']:
    deps[CFG[key]] = pa.file_hash(CFG[key])
pa.begin('terrazas')
terrazas_originales = pd.read_excel(CFG['terrazas_excel'], sheet_name='Hoja1')
revision_muestras = pd.read_excel(CFG['revision_muestras_excel'], sheet_name='Resumen')
seleccion_excel = pd.read_excel(CFG['seleccion_excel'], sheet_name='Hoja1').dropna(how='all')
chi_por_rio = {rid: pa.chi_1d(p) for rid,p in perfiles.items()}


# Notebook cell 5
def proyectar_terraza(p, este, norte):
    xy = p[["x_utm", "y_utm"]].to_numpy(float)
    inicio = xy[:-1]
    vector = xy[1:] - inicio
    longitud2 = np.sum(vector**2, axis=1)

    if np.any(longitud2 <= 0):
        raise ValueError("El cauce contiene segmentos espaciales de longitud cero.")

    punto = np.array([este, norte], dtype=float)
    fraccion = np.clip(
        np.sum((punto - inicio) * vector, axis=1) / longitud2,
        0, 1,
    )

    proyecciones = inicio + fraccion[:, None] * vector
    distancias = np.linalg.norm(proyecciones - punto, axis=1)
    j = int(np.argmin(distancias))
    f = fraccion[j]

    distancia = p["distancia_m"].to_numpy(float)
    elevacion = p["elevacion_m"].to_numpy(float)

    return {
        "distancia_proyectada_m": distancia[j] + f * (
            distancia[j + 1] - distancia[j]
        ),
        "distancia_lateral_m": distancias[j],
        "x_proyectado": proyecciones[j, 0],
        "y_proyectado": proyecciones[j, 1],
        "z_cauce_proyectado_m": elevacion[j] + f * (
            elevacion[j + 1] - elevacion[j]
        ),
    }


# Notebook cell 6
requeridas = {'Cauce', 'Cauce asociado', 'Terraza', 'Elevación (m)', 'Error elevación (m)',
    'Este UTM18S (m)', 'Norte UTM18S (m)', 'CRS coordenadas UTM'}
if requeridas - set(terrazas_originales): raise ValueError(f'Faltan columnas: {requeridas-set(terrazas_originales)}')
nombre_a_id = {pa.normalizar_nombre(v): k for k,v in nombres.items()}
if len(nombre_a_id) != len(nombres): raise ValueError('Nombres normalizados ambiguos.')
t = terrazas_originales.copy()
t['fila_excel'] = np.arange(len(t)) + 2
t['Cauce'] = t['Cauce'].ffill()
t['cauce_excel'] = t['Cauce asociado'].fillna(t['Cauce'])
t['ID'] = t.cauce_excel.map(lambda x: nombre_a_id.get(pa.normalizar_nombre(x)))
t['terraza_id'] = t.Terraza.astype('string').str.strip().str.replace('-', '_', regex=False)
for col in ['Elevación (m)', 'Error elevación (m)', 'Este UTM18S (m)', 'Norte UTM18S (m)']:
    t[col] = pd.to_numeric(t[col], errors='coerce')
filas, apartadas = [], []
for row in t.to_dict('records'):
    try:
        rid = row['ID']
        if rid not in perfiles or pd.isna(row['terraza_id']): raise ValueError('Sin asociación o identificador')
        if '32718' not in str(row['CRS coordenadas UTM']): raise ValueError('CRS de coordenadas no confirmado')
        if not np.isfinite([row[k] for k in ['Elevación (m)', 'Este UTM18S (m)', 'Norte UTM18S (m)']]).all():
            raise ValueError('Cota o coordenadas ausentes')
        row['asociacion_id'] = f"R{int(rid):02d}_{row['terraza_id']}"
        row.update(proyectar_terraza(perfiles[rid], row['Este UTM18S (m)'], row['Norte UTM18S (m)']))
        row['chi_proyectado_m'] = float(np.interp(row['distancia_proyectada_m'], perfiles[rid].distancia_m, chi_por_rio[rid]))
        row['proyeccion_en_extremo'] = bool(np.isclose(row['distancia_proyectada_m'], 0) or np.isclose(row['distancia_proyectada_m'], perfiles[rid].distancia_m.iloc[-1]))
        row['estado'] = 'ok'; filas.append(row)
    except (ValueError, KeyError, TypeError) as exc:
        row.update(estado='pendiente', detalle=str(exc)); apartadas.append(row)
terrazas = pd.DataFrame(filas, columns=[*t.columns, 'asociacion_id', 'distancia_proyectada_m',
    'distancia_lateral_m', 'x_proyectado', 'y_proyectado', 'z_cauce_proyectado_m',
    'chi_proyectado_m', 'proyeccion_en_extremo', 'estado'])
if terrazas.asociacion_id.duplicated().any(): raise ValueError('Asociación río-terraza duplicada; revisar la planilla.')
terrazas_no_asociadas = pd.DataFrame(apartadas, columns=[*t.columns, 'estado', 'detalle'])


# Notebook cell 8
filas_cruces, filas_asociaciones = [], []
for row in terrazas.to_dict('records'):
    rid, p = row['ID'], perfiles[row['ID']]
    cota, error = float(row['Elevación (m)']), row['Error elevación (m)']
    niveles = [('central', cota)]
    if pd.notna(error) and np.isfinite(error) and error >= 0:
        niveles += [('inferior', cota-error), ('superior', cota+error)]
    n_central = n_planos = 0
    for nivel, altura in niveles:
        posiciones, planos = pa.crossings(p, altura)
        for j, s in enumerate(posiciones, 1):
            filas_cruces.append({'ID': rid, 'terraza_id': row['terraza_id'], 'asociacion_id': row['asociacion_id'],
                'nivel': nivel, 'cruce_id': j, 'tipo': 'puntual', 'desde_m': s, 'hasta_m': s,
                'cota_m': altura, 'chi_m': float(np.interp(s, p.distancia_m, chi_por_rio[rid]))})
        for j, (a,b) in enumerate(planos, 1):
            filas_cruces.append({'ID': rid, 'terraza_id': row['terraza_id'], 'asociacion_id': row['asociacion_id'],
                'nivel': nivel, 'cruce_id': j, 'tipo': 'plano_ambiguo', 'desde_m': a, 'hasta_m': b,
                'cota_m': altura, 'chi_m': np.nan})
        if nivel == 'central': n_central, n_planos = len(posiciones), len(planos)
    row.update(n_cruces=n_central, n_planos=n_planos,
        estado_cruces='ambiguo' if n_planos or n_central>1 else 'unico' if n_central==1 else 'sin_cruce',
        perfil_hash=pa.profile_hash(p))
    filas_asociaciones.append(row)
asociaciones = pd.DataFrame(filas_asociaciones, columns=[*terrazas.columns, 'n_cruces', 'n_planos', 'estado_cruces', 'perfil_hash'])
cruces = pd.DataFrame(filas_cruces, columns=['ID','terraza_id','asociacion_id','nivel','cruce_id','tipo','desde_m','hasta_m','cota_m','chi_m'])


# Notebook cell 10
if {'id','Cauce','Terraza','ID knickpoint'} - set(seleccion_excel):
    raise ValueError('Faltan columnas de selección en el Excel.')
seleccion_excel['ID'] = seleccion_excel['id'].map(pa.normalize_id)
if seleccion_excel.ID.duplicated().any(): raise ValueError('Más de una selección por río en Excel.')
selecciones = []
for rid in pa.EXPECTED:
    ex = seleccion_excel.loc[seleccion_excel.ID.eq(rid)]
    kp = manuales.loc[manuales.ID.eq(rid)]
    fila = {'ID': rid, 'rio': nombres[rid], 'terraza_id': '', 'id_kp': '', 'id_kp_excel': '',
        'discrepancia_kp_excel': False, 'estado': 'pendiente', 'detalle': '', 'perfil_hash': pa.profile_hash(perfiles[rid])}
    motivos = []
    if len(kp) == 1: fila['id_kp'] = kp.iloc[0].id_kp
    else: motivos.append('sin_manual_valido' if kp.empty else 'multiples_manuales_requieren_eleccion')
    if ex.empty: motivos.append('sin_fila_excel')
    else:
        e = ex.iloc[0]
        if pa.normalizar_nombre(e.Cauce) != pa.normalizar_nombre(nombres[rid]): motivos.append('nombre_excel_incompatible')
        tid = pa.text(e.Terraza).replace('-', '_')
        if tid and not tid.lower().startswith('no tiene terraza'):
            fila['terraza_id'] = tid
            a = asociaciones.loc[asociaciones.ID.eq(rid) & asociaciones.terraza_id.eq(tid)]
            if len(a) != 1: motivos.append('terraza_sin_asociacion_unica')
        else: motivos.append('sin_terraza_seleccionada')
        fila['id_kp_excel'] = pa.text(e['ID knickpoint'])
        fila['discrepancia_kp_excel'] = bool(fila['id_kp_excel'] and fila['id_kp_excel'] != fila['id_kp'])
    fila.update(estado='pendiente' if motivos else 'seleccion_disponible', detalle='; '.join(motivos))
    selecciones.append(fila)
selecciones = pd.DataFrame(selecciones)
archivos = {k: pa.save_table('terrazas', k, table) for k,table in {
    'asociaciones': asociaciones, 'cruces': cruces, 'selecciones': selecciones,
    'terrazas_apartadas': terrazas_no_asociadas, 'manuales_rechazados': manuales_rechazados,
    'revision_muestras_fuente': revision_muestras}.items()}
pa.finish('terrazas', archivos, deps, extra={'criterio': 'Todas las asociaciones; ningún cruce elegido automáticamente'})
display(selecciones)


# Notebook cell 11
# ============================================================
# RECUPERAR LOS KP MANUALES PARA LOS GRÁFICOS
# ============================================================
# Lee tus selecciones originales y comprueba su ubicación
# contra los perfiles actuales.
# No modifica el CSV ni la tabla manuales usada en el análisis.

import numpy as np
import pandas as pd
from IPython.display import display

TOLERANCIA_XY_M = 0.10
TOLERANCIA_ELEVACION_M = 0.10

selecciones_originales = pd.read_csv(
    CFG["kp_manual_csv"],
    sep=CFG["kp_separador"],
    decimal=CFG["kp_decimal"],
    dtype={"ID": str},
).dropna(how="all")

selecciones_originales.columns = (
    selecciones_originales.columns.str.strip()
)

kp_recuperados = []
revision_kp = []

for fila in selecciones_originales.to_dict("records"):
    identificador = fila.get("id_kp")

    # Una fila sin identificador no es una selección de KP.
    if pd.isna(identificador) or not str(identificador).strip():
        continue

    rid = str(fila.get("ID", "")).strip()
    registro = {
        "ID": rid,
        "KP": str(identificador),
        "estado": "",
        "separacion_XY_m": np.nan,
        "diferencia_elevacion_m": np.nan,
    }

    try:
        rid = pa.normalize_id(fila["ID"])
        registro["ID"] = rid

        if rid not in perfiles:
            raise ValueError("Río no cargado")

        p = perfiles[rid]

        xy_original = np.array(
            [fila["x_utm"], fila["y_utm"]],
            dtype=float,
        )
        z_original = float(fila["elevacion_m"])

        if not np.isfinite([*xy_original, z_original]).all():
            raise ValueError("Selección sin coordenadas o elevación válidas")

        xy_actual = p[["x_utm", "y_utm"]].to_numpy(dtype=float)
        separaciones = np.linalg.norm(
            xy_actual - xy_original,
            axis=1,
        )

        # Busca la posición por coordenadas; no reutiliza
        # automáticamente el índice guardado anteriormente.
        candidatos = np.flatnonzero(
            separaciones <= TOLERANCIA_XY_M
        )

        if len(candidatos) != 1:
            raise ValueError(
                f"Se encontraron {len(candidatos)} nodos coincidentes; "
                "requiere revisión"
            )

        j = int(candidatos[0])
        nodo = p.iloc[j]
        diferencia_z = abs(
            float(nodo["elevacion_m"]) - z_original
        )

        registro["separacion_XY_m"] = float(separaciones[j])
        registro["diferencia_elevacion_m"] = diferencia_z

        if diferencia_z > TOLERANCIA_ELEVACION_M:
            raise ValueError("La elevación cambió; requiere revisión")

        recuperado = fila.copy()
        recuperado.update(
            ID=rid,
            indice=j,
            distancia_m=float(nodo["distancia_m"]),
            elevacion_m=float(nodo["elevacion_m"]),
            x_utm=float(nodo["x_utm"]),
            y_utm=float(nodo["y_utm"]),
        )

        kp_recuperados.append(recuperado)
        registro["estado"] = "Coincide con el perfil actual"

    except (ValueError, KeyError, TypeError) as exc:
        registro["estado"] = str(exc)

    revision_kp.append(registro)

manuales_grafico = pd.DataFrame(
    kp_recuperados,
    columns=selecciones_originales.columns,
)

print(f"KP disponibles para mostrar: {len(manuales_grafico)}")
display(pd.DataFrame(revision_kp))

# Notebook cell 13
# ============================================================
# PERFILES CON TERRAZAS Y KNICKPOINTS SELECCIONADOS MANUALMENTE
# ============================================================
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

# ---------- OPCIONES QUE PUEDES MODIFICAR ----------
RIOS_VER = []                  # [] = todos los ríos cargados
TIPO_PERFIL = "chi"               # "chi" o "longitudinal"

MOSTRAR_TERRAZAS = True
MOSTRAR_KP = True

ETIQUETAR_TERRAZAS = True         # Nombre junto a cada terraza
ETIQUETAR_KP = True               # Identificador junto a cada KP

GUARDAR_FIGURAS = False
# --------------------------------------------------

tipo = TIPO_PERFIL.strip().lower()
if tipo not in {"chi", "longitudinal"}:
    raise ValueError('TIPO_PERFIL debe ser "chi" o "longitudinal".')

rios_elegidos = (
    list(map(str, RIOS_VER))
    if RIOS_VER
    else list(map(str, perfiles.keys()))
)

# La fecha corresponde al momento en que ejecutas esta celda.
fecha_figuras = datetime.now().strftime("%Y%m%d_%H%M%S_%f")

for rid in rios_elegidos:
    if rid not in perfiles:
        print(f"Río {rid}: no está cargado en perfiles.")
        continue

    p = perfiles[rid]
    distancia = p["distancia_m"].to_numpy(dtype=float)
    elevacion = p["elevacion_m"].to_numpy(dtype=float)
    nombre = nombres.get(rid, f"Río {rid}")

    if tipo == "chi":
        coordenada_x = np.asarray(chi_por_rio[rid], dtype=float)
        etiqueta_x = "χ (m)"
        titulo_tipo = "Perfil χ–elevación"
    else:
        coordenada_x = distancia / 1000.0
        etiqueta_x = "Distancia desde el outlet recortado (km)"
        titulo_tipo = "Perfil longitudinal"

    fig, ax = plt.subplots(figsize=(11, 5))

    ax.plot(
        coordenada_x,
        elevacion,
        color="steelblue",
        linewidth=1.8,
        label="Perfil del río",
        zorder=2,
    )

    # ---------- TERRAZAS ----------
    # Se muestra la cota de la terraza en su posición proyectada
    # sobre el río. Esta posición NO es un cruce por cota.
    if MOSTRAR_TERRAZAS:
        terrazas_rio = asociaciones.loc[
            asociaciones["ID"].astype(str).eq(rid)
        ]

        for t in terrazas_rio.to_dict("records"):
            if tipo == "chi":
                x_terraza = float(t["chi_proyectado_m"])
            else:
                x_terraza = float(t["distancia_proyectada_m"]) / 1000.0

            z_terraza = float(t["Elevación (m)"])
            nombre_terraza = str(t["terraza_id"])

            if not np.isfinite([x_terraza, z_terraza]).all():
                continue

            ax.scatter(
                x_terraza,
                z_terraza,
                marker="^",
                s=75,
                edgecolors="black",
                linewidths=0.5,
                label=f"Terraza {nombre_terraza}",
                zorder=4,
            )

            if ETIQUETAR_TERRAZAS:
                ax.annotate(
                    nombre_terraza,
                    (x_terraza, z_terraza),
                    xytext=(5, 7),
                    textcoords="offset points",
                    fontsize=8,
                )

    # ---------- KP MANUALES ----------
    # Usa únicamente la tabla de selecciones manuales válidas
    # cargada anteriormente mediante pa.load_manuals(perfiles).
    if MOSTRAR_KP:
        kp_rio = manuales_grafico.loc[manuales_grafico["ID"].astype(str).eq(rid)]

        if kp_rio.empty:
            print(f"{rid} · {nombre}: sin KP manual válido cargado.")

        for numero, k in enumerate(kp_rio.to_dict("records")):
            indice = int(k["indice"])

            if not 0 <= indice < len(p):
                raise ValueError(
                    f"KP {k['id_kp']}: índice fuera del perfil "
                    f"del río {rid}. Revisa la selección."
                )

            x_kp = coordenada_x[indice]
            z_kp = elevacion[indice]

            ax.scatter(
                x_kp,
                z_kp,
                marker="*",
                color="black",
                edgecolors="white",
                linewidths=0.7,
                s=190,
                label="KP seleccionado manualmente" if numero == 0 else None,
                zorder=6,
            )

            if ETIQUETAR_KP:
                ax.annotate(
                    str(k["id_kp"]),
                    (x_kp, z_kp),
                    xytext=(7, -15),
                    textcoords="offset points",
                    fontsize=8,
                    color="black",
                )

    ax.set(
        title=f"{rid} · {nombre} — {titulo_tipo}",
        xlabel=etiqueta_x,
        ylabel="Elevación (m)",
    )
    ax.grid(alpha=0.2)

    # Mantiene la leyenda de los KP incluso sin mostrar terrazas.
    handles, labels = ax.get_legend_handles_labels()
    leyenda = dict(zip(labels, handles))
    ax.legend(
        leyenda.values(),
        leyenda.keys(),
        fontsize=8,
        loc="best",
    )

    fig.tight_layout()

    if GUARDAR_FIGURAS:
        archivo = pa.output(
            "terrazas",
            f"fig_rio_{int(rid):02d}_{tipo}_{fecha_figuras}",
            "png",
        )
        fig.savefig(archivo, dpi=300, bbox_inches="tight")
        print(f"Figura guardada: {archivo}")

    plt.show()
    plt.close(fig)


# Notebook cell 14

