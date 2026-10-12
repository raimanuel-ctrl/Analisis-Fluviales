# Code-only reading snapshot; not an independently executable script.
# Notebook is authoritative. Cell magics, external data and cell order still apply.

# Notebook cell 3
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json, hashlib
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from IPython.display import display

BASE = Path(r'C:\datos\Memoria\Codigos Stanford')
RUN = '20261005_post_auditoria_01'
MANIFEST = BASE / 'Productos' / 'post_auditoria' / RUN / 'perfiles' / f'manifest_{RUN}.json'
FUENTE_MUESTRAS = BASE / 'Productos' / 'post_auditoria' / RUN / 'terrazas' / 'revision_muestras_fuente_20261005_post_auditoria_01_20261005_152125_119188.csv'
FUENTE_TERRAZAS = Path(r'C:\datos\Memoria_Codex\outputs\terrazas_audit_20260923\Terrazas_auditadas.csv')
THETA, A0, N = 0.40, 1e6, 1.0
MIN_NODOS = 5  # Se cuentan puntos originales; reducirlo aumenta la fragilidad del ajuste.
GUARDAR_TABLAS = True
GUARDAR_FIGURAS = False
FECHA = datetime.now(ZoneInfo('America/Los_Angeles')).strftime('%Y%m%d_%H%M%S_%f')
SALIDA = BASE / 'Productos' / 'K_seleccion' / FECHA

MUESTRAS = {
    'CONT 11': dict(edad_ka=124.1, error_ka=10.8, edad_asignada_ka=125., OD=11.70, U_publicada_m_ka=0.54),
    'PUHE 11': dict(edad_ka=111.8, error_ka=9.5, edad_asignada_ka=125., OD=8.98, U_publicada_m_ka=0.53),
    'PUHE 10': dict(edad_ka=102.3, error_ka=9.3, edad_asignada_ka=106., OD=14.8, U_publicada_m_ka=0.50),
}
RIOS = {
    '3': dict(nombre='Estero Quebrada Honda', muestra='CONT 11', terraza='JM2015_028', cota=67.80),
    '4': dict(nombre='Estero Malhueco', muestra='CONT 11', terraza='JM2015_030', cota=85.60),
    '7': dict(nombre='Estero de Chanco', muestra='PUHE 11', terraza='JM2015_045', cota=62.68),
    '8': dict(nombre='Río Rahue', muestra='PUHE 11', terraza='JM2015_046', cota=93.47),
    '9': dict(nombre='Río Curanilahue', muestra='PUHE 11', terraza='JM2015_048', cota=70.74),
    '10': dict(nombre='Río Curanipe', muestra='PUHE 10', terraza='JM2015_163', cota=24.50),
}
RIOS_VER = ['3', '4', '7', '8', '9', '10']
CRUCES_ELEGIDOS = {}  # Ejemplo: {'8': 2}. Vacío: evaluar todos los cruces como alternativas.
ALUVIALES = {rid: None for rid in RIOS}
EVIDENCIA_ALUVIAL = {rid: '' for rid in RIOS}
# Ejemplo de formato (NO son límites reales):
# ALUVIALES['8'] = [(0., 300.), (700., 850.)]
# EVIDENCIA_ALUVIAL['8'] = 'Cartografía/unidad, observación y fecha que sustentan estos límites'
# ALUVIALES['9'] = []  # Solo si se revisó y no hay un tramo que excluir.
# EVIDENCIA_ALUVIAL['9'] = 'Revisado con ...; no se identificó planicie aluvial en el intervalo'

# Corrección vertical EXPLÍCITA del escenario de relieve/edad, respecto al outlet actual.
# Cero es una hipótesis exploratoria, NO un nivel marino antiguo demostrado.
CORRECCION_VERTICAL_M = {rid: 0. for rid in RIOS}
display(pd.DataFrame.from_dict(RIOS, orient='index').rename_axis('ID'))
print('Selección cargada: seis ríos; Chanco se conserva como prueba.')


# Notebook cell 5
manifest_bytes = MANIFEST.read_bytes()
manifest = json.loads(manifest_bytes.decode('utf-8-sig'))
fuentes_documentales = []
for archivo_fuente in [FUENTE_MUESTRAS, FUENTE_TERRAZAS]:
    fuentes_documentales.append(dict(path=str(archivo_fuente),
        sha256=hashlib.sha256(archivo_fuente.read_bytes()).hexdigest() if archivo_fuente.exists() else None))
if manifest.get('estado') != 'completo':
    raise ValueError('La exportación de perfiles no está completa.')
perfiles_K, qc, entradas = {}, [], []
for rid, cfg in RIOS.items():
    try:
        item = manifest['files'][f'rio_{rid}']
        archivo = Path(item['path'])
        huella = hashlib.sha256(archivo.read_bytes()).hexdigest()
        if huella != item['sha256']:
            raise ValueError('El CSV cambió después de su exportación.')
        p = pd.read_csv(archivo, float_precision='round_trip')
        columnas = ['distancia_m', 'elevacion_m', 'area_m2', 'x_utm', 'y_utm']
        if set(columnas) - set(p):
            raise ValueError('Faltan columnas de geometría del perfil.')
        s, z, a = (p[c].to_numpy(float) for c in columnas[:3])
        if len(s) < 2 or not np.isfinite(s).all() or not np.isclose(s[0], 0.) or np.any(np.diff(s) <= 0):
            raise ValueError('Distancias inválidas: se esperaba cero en el outlet y orden creciente.')
        if not np.isfinite(a).all() or np.any(a <= 0) or not np.isfinite(z[0]):
            raise ValueError('Área o elevación del outlet inválida.')
        p['chi_040_m'] = cumulative_trapezoid((A0 / a)**THETA, x=s, initial=0.)
        perfiles_K[rid] = p
        entradas.append(dict(ID=rid, path=str(archivo), sha256=huella))
        qc.append(dict(ID=rid, rio=cfg['nombre'], estado='cargado', nodos=len(p),
                       z_outlet_m=z[0], huecos_z=int((~np.isfinite(z)).sum()), detalle=''))
    except (OSError, KeyError, ValueError) as exc:
        qc.append(dict(ID=rid, rio=cfg['nombre'], estado='error', detalle=str(exc)))
qc_K = pd.DataFrame(qc)
display(qc_K)
if not perfiles_K:
    raise RuntimeError('No se cargó ningún perfil; revisar las rutas y los errores anteriores.')


# Notebook cell 7
def cruces_cota(p, cota):
    s, z = p.distancia_m.to_numpy(), p.elevacion_m.to_numpy()
    puntos, planos = [], []
    for j in range(len(s)-1):
        if not np.isfinite(z[j:j+2]).all():
            continue
        if z[j] == z[j+1]:
            if z[j] == cota:
                planos.append([float(s[j]), float(s[j+1])])
        elif min(z[j], z[j+1]) <= cota <= max(z[j], z[j+1]):
            puntos.append(float(s[j] + (cota-z[j])*(s[j+1]-s[j])/(z[j+1]-z[j])))
    # Agrupar duplicados en vértices sin redondear límites fuera del dominio.
    unicos = []
    for v in sorted(puntos):
        if not unicos or abs(v-unicos[-1]) > 1e-6:
            unicos.append(v)
    return unicos, planos

filas, limites = [], {}
for rid, p in perfiles_K.items():
    puntos, planos = cruces_cota(p, RIOS[rid]['cota'])
    limites[rid] = (puntos, planos)
    if not puntos or planos:
        filas.append(dict(ID=rid, cruce_id=None, distancia_m=np.nan,
                          estado='plano_a_cota' if planos else 'sin_cruce', planos=str(planos)))
    for i, s in enumerate(puntos, 1):
        filas.append(dict(ID=rid, cruce_id=i, distancia_m=s,
                          estado='alternativa' if len(puntos)>1 else 'unico', planos=str(planos)))
cruces_K = pd.DataFrame(filas, columns=['ID','cruce_id','distancia_m','estado','planos'])
display(cruces_K)


# Notebook cell 9
def intervalos_roca(fin, cortes):
    normalizados = []
    for a, b in cortes:
        a, b = float(a), float(b)
        if not np.isfinite([a,b]).all() or a < 0 or b <= a:
            raise ValueError('Intervalo aluvial inválido; usar 0 <= inicio < fin.')
        if a < fin:
            normalizados.append((a, min(b, fin)))
    unidos = []
    for a,b in sorted(normalizados):
        if unidos and a <= unidos[-1][1]:
            unidos[-1] = (unidos[-1][0], max(b, unidos[-1][1]))
        else:
            unidos.append((a,b))
    tramos, actual = [], 0.
    for a,b in unidos:
        if a > actual:
            tramos.append((actual,a))
        actual = b
    if actual < fin:
        tramos.append((actual,fin))
    if not tramos:
        raise ValueError('Las exclusiones eliminan todo el intervalo.')
    return tramos, unidos

def ajuste_tramos(p, tramos):
    s, z, chi = (p[c].to_numpy(float) for c in ['distancia_m','elevacion_m','chi_040_m'])
    grupos, detalles = [], []
    for a,b in tramos:
        if not 0 <= a < b <= s[-1]:
            raise ValueError('Tramo fuera del dominio o de longitud cero.')
        lo = max(0, np.searchsorted(s,a,side='right')-1)
        hi = min(len(s)-1,np.searchsorted(s,b))
        if not np.isfinite(z[lo:hi+1]).all():
            raise ValueError('Hueco de elevación dentro del tramo; no se interpola sobre él.')
        dentro = (s>=a)&(s<=b)
        if dentro.sum() < MIN_NODOS:
            raise ValueError(f'Tramo {a:.1f}–{b:.1f} m: menos de {MIN_NODOS} nodos originales.')
        q = np.unique(np.r_[a,s[dentro],b])
        x, y = np.interp(q,s,chi), np.interp(q,s,z)
        xc, yc = x-x.mean(), y-y.mean()
        if xc@xc <= 0:
            raise ValueError('Tramo sin rango de chi.')
        grupos.append((xc,yc))
        detalles.append(dict(inicio_m=a,fin_m=b,n_nodos=int(dentro.sum()),
                             Mchi_individual=float((xc@yc)/(xc@xc))))
    xx = sum(x@x for x,y in grupos)
    pendiente = sum(x@y for x,y in grupos)/xx
    if not np.isfinite(pendiente) or pendiente <= 0:
        raise ValueError('Pendiente común no positiva; no se usa valor absoluto.')
    ss = sum(y@y for x,y in grupos)
    r2 = 1-sum(np.sum((y-pendiente*x)**2) for x,y in grupos)/ss if ss>0 else np.nan
    return float(pendiente), float(r2), detalles

resultados, detalles_tramos = [], []
for rid, cfg in RIOS.items():
    m = MUESTRAS[cfg['muestra']]
    base = dict(ID=rid,rio=cfg['nombre'],prueba=(rid=='7'),muestra=cfg['muestra'],
                terraza=cfg['terraza'],cota_m=cfg['cota'],edad_ka=m['edad_ka'],
                error_edad_ka=m['error_ka'],OD_pct=m['OD'],theta=THETA,n=N,
                U_publicada_m_ka=m['U_publicada_m_ka'],correccion_m=CORRECCION_VERTICAL_M[rid],
                unidad_K='m^0.2/año',asociacion='correlacion_de_sector',K_relieve_edad=np.nan,K_U_publicada=np.nan)
    p = perfiles_K.get(rid)
    puntos, planos = limites.get(rid,([],[]))
    opciones = list(enumerate(puntos,1))
    if rid in CRUCES_ELEGIDOS:
        opciones = [(i,s) for i,s in opciones if i==CRUCES_ELEGIDOS[rid]]
    if p is None or planos or not opciones:
        for escenario in ['completo','sin_aluvion']:
            resultados.append(dict(base,escenario=escenario,estado='pendiente_geometria',
                                    detalle='Ver QC/cruces; falta perfil, hay segmento plano o cruce inexistente.'))
        continue
    for i,fin in opciones:
        for escenario in ['completo','sin_aluvion']:
            row = dict(base,cruce_id=i,fin_m=fin,n_cruces=len(puntos),escenario=escenario,
                       estado='pendiente',evidencia_aluvial=EVIDENCIA_ALUVIAL[rid])
            try:
                if not 0 <= m['OD'] < 15:
                    raise ValueError('Muestra fuera del filtro OD < 15%.')
                if escenario=='sin_aluvion' and (ALUVIALES[rid] is None or not EVIDENCIA_ALUVIAL[rid].strip()):
                    row.update(estado='pendiente_aluvial',detalle='Completar intervalos y evidencia; el escenario completo sí se evalúa.')
                    resultados.append(row)
                    continue
                cortes = [] if escenario=='completo' else ALUVIALES[rid]
                if any(float(b)>float(p.distancia_m.iloc[-1]) for a,b in cortes):
                    raise ValueError('Un corte aluvial queda fuera del perfil completo.')
                tramos, excluidos = intervalos_roca(fin,cortes)
                M, r2, det = ajuste_tramos(p,tramos)
                ksn = M*A0**THETA
                edad, error = m['edad_ka'],m['error_ka']
                if not np.isfinite([edad,error]).all() or edad<=error or error<0:
                    raise ValueError('Edad/error inválidos.')
                relieve = cfg['cota']-float(p.elevacion_m.iloc[0])-CORRECCION_VERTICAL_M[rid]
                if not np.isfinite(relieve) or relieve<=0:
                    raise ValueError('Relieve corregido no positivo; revisar datum y corrección.')
                row.update(estado='escenario_condicional',detalle='Alternativas de cruce' if len(puntos)>1 else 'Cruce único',
                           Mchi=M,ksn_040=ksn,R2_dentro_tramos=r2,n_tramos=len(tramos),
                           longitud_roca_m=sum(b-a for a,b in tramos),
                           longitud_excluida_m=sum(b-a for a,b in excluidos),
                           cortes_aplicados=json.dumps(excluidos),relieve_corregido_m=relieve,
                           E_aparente_m_ka=relieve/edad,
                           K_relieve_edad=relieve/(edad*1000*ksn),
                           K_edad_mas_error=relieve/((edad+error)*1000*ksn),
                           K_edad_menos_error=relieve/((edad-error)*1000*ksn),
                           K_U_publicada=m['U_publicada_m_ka']/(1000*ksn))
                for d in det:
                    detalles_tramos.append(dict(ID=rid,cruce_id=i,escenario=escenario,**d))
            except (ValueError,TypeError) as exc:
                row.update(estado='no_calculable',detalle=str(exc))
            resultados.append(row)
resultados_K = pd.DataFrame(resultados)
tramos_K = pd.DataFrame(detalles_tramos)
display(resultados_K)
display(tramos_K)
print('Estado de cada escenario:')
display(resultados_K.groupby(['escenario','estado']).size().rename('filas').reset_index())


# Notebook cell 11
if GUARDAR_TABLAS:
    SALIDA.mkdir(parents=True,exist_ok=True)
    for nombre,tabla in [('resultados_K',resultados_K),('tramos_K',tramos_K),('cruces_K',cruces_K),('qc_K',qc_K)]:
        tabla.to_csv(SALIDA/f'{nombre}_{FECHA}.csv',index=False,encoding='utf-8-sig')
    for rid,p in perfiles_K.items():
        p[['distancia_m','elevacion_m','area_m2','x_utm','y_utm','chi_040_m']].to_csv(
            SALIDA/f'perfil_theta040_{rid}_{FECHA}.csv',index=False,encoding='utf-8-sig')
    registro=dict(fecha=FECHA,zona='America/Los_Angeles',theta=THETA,A0_m2=A0,n=N,
                  manifest=str(MANIFEST),manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                  fuentes_documentales=fuentes_documentales,
                  entradas=entradas,rios=RIOS,muestras=MUESTRAS,cruces_elegidos=CRUCES_ELEGIDOS,
                  aluviales=ALUVIALES,evidencia=EVIDENCIA_ALUVIAL,correccion_vertical_m=CORRECCION_VERTICAL_M,
                  min_nodos=MIN_NODOS,metodo='pendiente comun, interceptos independientes por tramo',
                  alcance='K exploratorio condicionado; equilibrio y correlaciones no demostrados')
    (SALIDA/f'parametros_K_{FECHA}.json').write_text(json.dumps(registro,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Archivos guardados en:',SALIDA)
else:
    print('Tablas disponibles en memoria; GUARDAR_TABLAS=False.')


# Notebook cell 13
get_ipython().run_line_magic('matplotlib','inline')
import matplotlib.pyplot as plt
GROSOR = 1.8
BORDE_NEGRO = False
RELLENO = False
GRILLA = True
MAX_PUNTOS_GRAFICO = 2500
# BORDE_NEGRO = True
# GROSOR = 2.5
# RELLENO = True
# GRILLA = False

for rid in RIOS_VER:
    if rid not in perfiles_K:
        continue
    p,cfg = perfiles_K[rid],RIOS[rid]
    indices = np.unique(np.linspace(0,len(p)-1,min(len(p),MAX_PUNTOS_GRAFICO)).astype(int))
    fig,ax = plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
    for eje,col in zip(ax,['distancia_m','chi_040_m']):
        x=p[col].to_numpy()/1000
        z=p.elevacion_m.to_numpy()
        if BORDE_NEGRO:
            eje.plot(x[indices],z[indices],color='black',lw=GROSOR+1.5)
        eje.plot(x[indices],z[indices],color='#24788E',lw=GROSOR,label='Perfil')
        if RELLENO:
            eje.fill_between(x[indices],z[indices],float(z[0]),color='#24788E',alpha=.15)
        eje.axhline(cfg['cota'],color='#AD5826',lw=1.2,label=cfg['terraza'])
        for i,s in enumerate(limites[rid][0],1):
            v=np.interp(s,p.distancia_m,p[col])/1000
            eje.scatter([v],[cfg['cota']],color='#AD5826',zorder=5)
            eje.annotate(str(i),(v,cfg['cota']),xytext=(5,5),textcoords='offset points')
        for a,b in ALUVIALES[rid] or []:
            eje.axvspan(np.interp(a,p.distancia_m,p[col])/1000,
                        np.interp(b,p.distancia_m,p[col])/1000,color='gray',alpha=.2)
        eje.set(xlabel='Distancia desde outlet (km)' if col=='distancia_m' else 'χ (km)',ylabel='Elevación (m)')
        eje.grid(GRILLA,alpha=.2)
    ax[0].legend(fontsize=8)
    fig.suptitle(f"{rid} · {cfg['nombre']} · {cfg['muestra']}" + (' · PRUEBA' if rid=='7' else ''))
    if GUARDAR_FIGURAS:
        SALIDA.mkdir(parents=True,exist_ok=True)
        fig.savefig(SALIDA/f'perfil_{rid}_{FECHA}.png',dpi=200)
    plt.show()
    plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
for rid in RIOS_VER:
    if rid not in perfiles_K:
        continue
    p=perfiles_K[rid]
    indices=np.unique(np.linspace(0,len(p)-1,min(len(p),MAX_PUNTOS_GRAFICO)).astype(int))
    z=p.elevacion_m.to_numpy()-float(p.elevacion_m.iloc[0])
    for eje,col in zip(ax,['distancia_m','chi_040_m']):
        eje.plot(p[col].to_numpy()[indices]/1000,z[indices],lw=1.2,label=RIOS[rid]['nombre'])
for eje,label in zip(ax,['Distancia desde outlet (km)','χ desde outlet (km)']):
    eje.set(xlabel=label,ylabel='Elevación sobre outlet (m)',xlim=(0,None))
    eje.grid(GRILLA,alpha=.2)
ax[1].legend(fontsize=7)
if GUARDAR_FIGURAS:
    SALIDA.mkdir(parents=True,exist_ok=True)
    fig.savefig(SALIDA/f'perfiles_origen_comun_{FECHA}.png',dpi=200)
plt.show()
plt.close(fig)

