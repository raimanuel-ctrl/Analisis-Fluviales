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

pa.NOTEBOOK_PATH = CARPETA_NOTEBOOKS / 'Interpretaciones_K_deltaT_20261005.ipynb'


# Notebook cell 3
perfiles, nombres, manifest = pa.load_profiles()
selecciones = pa.stage_table('terrazas', 'selecciones')
asociaciones = pa.stage_table('terrazas', 'asociaciones')
cruces = pa.stage_table('terrazas', 'cruces')
manuales, manuales_rechazados = pa.load_manuals(perfiles)
deps = {**pa.dependency('perfiles'), **pa.dependency('terrazas')}
deps[CFG['kp_manual_csv']] = pa.file_hash(CFG['kp_manual_csv'])
pa.begin('interpretaciones')


# Notebook cell 5
columnas = ['ID','escenario','terraza_id','id_kp','cruce_id','inicio_roca_m','U_m_ka','fuente_U',
    'muestra_respaldo','edad_adoptada_ka','fuente_edad','justificacion_asociacion',
    'aplicar_equilibrio','justificacion_equilibrio','transferir_K','justificacion_transferencia',
    'geologia_revisada','theta','A0_m2','n']
plantilla = pd.DataFrame([{'ID': rid, 'escenario': 'base', 'terraza_id': '', 'id_kp': '',
    'cruce_id': np.nan, 'inicio_roca_m': np.nan, 'U_m_ka': np.nan, 'fuente_U': '',
    'muestra_respaldo': '', 'edad_adoptada_ka': np.nan, 'fuente_edad': '', 'justificacion_asociacion': '',
    'aplicar_equilibrio': False, 'justificacion_equilibrio': '', 'transferir_K': False,
    'justificacion_transferencia': '', 'geologia_revisada': False,
    'theta': CFG['theta'], 'A0_m2': CFG['A0_m2'], 'n': 1} for rid in pa.EXPECTED], columns=columnas)
ruta_plantilla = pa.save_table('interpretaciones', 'plantilla_escenarios', plantilla)
if CFG['escenarios_csv']:
    escenarios = pd.read_csv(CFG['escenarios_csv'], dtype={'ID': str})
    if set(columnas)-set(escenarios): raise ValueError(f'Faltan columnas: {set(columnas)-set(escenarios)}')
    escenarios['ID'] = escenarios.ID.map(pa.normalize_id)
    if escenarios[['ID','escenario']].duplicated().any(): raise ValueError('Escenarios duplicados.')
    if escenarios.escenario.isna().any(): raise ValueError('Faltan nombres de escenario.')
    faltan = set(pa.EXPECTED)-set(escenarios.ID)
    escenarios = pd.concat([escenarios, plantilla.loc[plantilla.ID.isin(faltan)]], ignore_index=True)
    deps[CFG['escenarios_csv']] = pa.file_hash(CFG['escenarios_csv'])
else:
    escenarios = plantilla.copy()
print('Plantilla de decisiones:', ruta_plantilla)


# Notebook cell 7
def afirmativo(value):
    return str(value).strip().lower() in {'true','1','1.0','si','sí'}

def evaluar_escenario(e):
    rid = e['ID']; p = perfiles[rid]
    r = {**e, 'run_id': RUN_ID, 'perfil_hash': pa.profile_hash(p), 'rio': nombres[rid],
         'estado': 'pendiente', 'detalle': '', 'ksn_calibracion': np.nan, 'K': np.nan,
         'unidad_K': '', 'delta_t_ka': np.nan, 'R2_calibracion': np.nan,
         'limite_inferior_m': np.nan, 'limite_superior_m': np.nan, 'posicion_kp_m': np.nan}
    try:
        seleccion = selecciones.loc[selecciones.ID.eq(rid)]
        if len(seleccion) != 1: raise ValueError('Selección regional ausente o duplicada')
        sel = seleccion.iloc[0]
        if 'nombre_excel_incompatible' in pa.text(sel.detalle): raise ValueError('Nombre del Excel incompatible')
        tid = pa.text(e['terraza_id'])
        if not tid or tid != pa.text(sel.terraza_id): raise ValueError('Escenario sin terraza vigente del Excel')
        asociado = asociaciones.loc[asociaciones.ID.eq(rid) & asociaciones.terraza_id.eq(tid)]
        if len(asociado) != 1: raise ValueError('Asociación río-terraza no única')
        kp = manuales.loc[manuales.ID.eq(rid) & manuales.id_kp.eq(pa.text(e['id_kp']))]
        if len(kp) != 1: raise ValueError('KP no elegido de forma única desde el CSV manual válido')
        j = float(e['cruce_id'])
        if not np.isfinite(j) or not j.is_integer(): raise ValueError('Falta cruce entero adoptado')
        opciones = cruces.loc[cruces.ID.eq(rid) & cruces.terraza_id.eq(tid) & cruces.nivel.eq('central')]
        if opciones.tipo.eq('plano_ambiguo').any(): raise ValueError('Cota sobre segmento plano: requiere revisión específica')
        c = opciones.loc[opciones.tipo.eq('puntual') & opciones.cruce_id.eq(int(j))]
        if len(c) != 1: raise ValueError('Cruce adoptado no disponible')
        inicio, corte, fin = float(e['inicio_roca_m']), float(c.iloc[0].desde_m), float(kp.iloc[0].distancia_m)
        if not np.isfinite([inicio,corte,fin]).all() or not 0 <= inicio < corte < fin <= p.distancia_m.iloc[-1]:
            raise ValueError('Se requiere inicio_roca < cruce < KP dentro del perfil')
        theta, A0, n = float(e['theta']), float(e['A0_m2']), float(e['n'])
        if n != 1 or not np.isclose(theta, CFG['theta']) or not np.isclose(A0, CFG['A0_m2']):
            raise ValueError('Este escenario requiere n=1 y theta/A0 de la corrida común')
        muestra = re.sub(r'\s+', ' ', pa.text(e['muestra_respaldo']).upper().replace('_',' '))
        if muestra not in CFG['muestras_autorizadas']: raise ValueError('Muestra no autorizada en la configuración')
        respaldo = re.sub(r'\s+', ' ', pa.text(asociado.iloc[0].get('Muestra(s) de respaldo de la correlación', '')).upper().replace('_',' '))
        if not re.search(r'(?<![A-Z0-9])'+re.escape(muestra)+r'(?![A-Z0-9])', respaldo):
            raise ValueError('La muestra no figura en el respaldo documentado de esta asociación')
        edad, U = float(e['edad_adoptada_ka']), float(e['U_m_ka'])
        if not np.isfinite([edad,U]).all() or edad<=0 or U<=0: raise ValueError('Faltan edad y U positivas con unidades')
        for field in ['fuente_U','fuente_edad','justificacion_asociacion','justificacion_equilibrio']:
            if not pa.text(e[field]): raise ValueError(f'Falta {field}')
        if not afirmativo(e['aplicar_equilibrio']): raise ValueError('Equilibrio no adoptado')
        ajuste = pa.interval_fit(p, inicio, corte, theta, A0)
        K = (U/1000)/ajuste['ksn']
        r.update(estado='K_condicional', detalle='Transferencia de K no adoptada',
            K=K, ksn_calibracion=ajuste['ksn'], R2_calibracion=ajuste['R2'],
            unidad_K=f'm^{1-2*theta:g}/año', limite_inferior_m=inicio,
            limite_superior_m=corte, posicion_kp_m=fin)
        if afirmativo(e['transferir_K']):
            if not afirmativo(e['geologia_revisada']) or not pa.text(e['justificacion_transferencia']):
                r['detalle'] = 'K calculado; transferencia pendiente de revisión geológica y justificación'
            else:
                chi = ajuste['chi']
                delta_chi = float(np.interp(fin,p.distancia_m,chi)-np.interp(corte,p.distancia_m,chi))
                r.update(delta_t_ka=delta_chi/(K*A0**theta)/1000, estado='K_y_tiempo_modelado',
                    detalle='n=1, K transferido constante, área estacionaria; no es edad tectónica')
    except (ValueError, TypeError, KeyError, IndexError) as exc:
        r['detalle'] = str(exc)
    return r

import re
resultados = pd.DataFrame([evaluar_escenario(e) for e in escenarios.to_dict('records')])
if set(resultados.ID) != set(pa.EXPECTED): raise ValueError('Faltan ríos en la evaluación regional.')
archivos = {'resultados': pa.save_table('interpretaciones','escenarios_K_deltaT',resultados),
    'plantilla': ruta_plantilla, 'escenarios_entrada': pa.save_table('interpretaciones','escenarios_entrada',escenarios)}
pa.finish('interpretaciones', archivos, deps, extra={'advertencia': 'Resultados condicionales; sin inferencia automática de edades tectónicas'})
display(resultados[['ID','rio','escenario','estado','K','delta_t_ka','detalle']])

