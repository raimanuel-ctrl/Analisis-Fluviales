"""Edades de highstand asignadas en la Tabla 2 publicada; política 20/10/2026."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
import numpy as np

URL='https://earthweb.ess.washington.edu/gomberg/SZOSeminar/Papers/MunozQuatSci15.pdf'
DOI='10.1016/j.quascirev.2015.01.005'
POLITICA={'fecha_decision':'2026-10-20','fuente':'Jara-Muñoz et al. (2015), Tabla 2, p. 179',
          'doi':DOI,'url':URL,'columna':'Assigned age (ka), no Age (ka) de la muestra',
          'edades_ka':{'5e':125.0,'5c':106.0},'rango_ka':{'5e':5.0,'5c':5.0},
          'datos_internos':'Geometría y correlación espacial; no determinan la edad numérica adoptada.',
          'alcance_rango':'Rango adoptado para highstand; no ocho dataciones independientes.'}

def aplicar(datos):
    d=datos.copy()
    unknown=set(d.mis.dropna().astype(str))-set(POLITICA['edades_ka'])
    if unknown:raise ValueError('MIS sin política cronológica publicada: '+str(unknown))
    if 'T_ka_original_catalogo' not in d:d['T_ka_original_catalogo']=d.T_ka
    if 'error_T_ka_original_catalogo' not in d:d['error_T_ka_original_catalogo']=d.error_T_ka
    d['T_ka']=d.mis.map(POLITICA['edades_ka']).astype(float)
    d['error_T_ka']=d.mis.map(POLITICA['rango_ka']).astype(float)
    d['fuente_edad_calibracion']=POLITICA['fuente']
    d['doi_edad_calibracion']=DOI
    d['tipo_edad_calibracion']='Edad de highstand asignada al marcador por correlación MIS; no edad directa de muestra o cauce.'
    d['edad_usada']=d['tipo_edad_calibracion']
    d['politica_cronologica']='JM2015_Tabla2_20261020'
    return d

def actualizar_resultados(datos):
    """Recalcular columnas de K derivadas de T, conservando posiciones y roles."""
    old=datos.copy();d=aplicar(datos)
    scale=old.T_ka.to_numpy(dtype=float)/d.T_ka.to_numpy(dtype=float)
    d['K_m02_anio']=(d.chi_PT_m-d.chi_B_m)/(1e6**.4*d.T_ka*1000)
    if 'direccion' in d:d.loc[d.direccion.eq('meseta'),'K_m02_anio']=np.nan
    for col in ['K_intervalo_min','K_intervalo_max','error_K_H_aprox']:
        if col in d:d[col]=old[col]*scale
    if 'error_K_T_aprox' in d:d['error_K_T_aprox']=d.K_m02_anio*d.error_T_ka/d.T_ka
    if 'error_K_propagado_aprox' in d:
        d['error_K_propagado_aprox']=np.hypot(d.error_K_H_aprox,d.error_K_T_aprox)
    return d
