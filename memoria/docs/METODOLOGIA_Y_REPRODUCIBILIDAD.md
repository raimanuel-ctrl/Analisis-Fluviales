# Qué reproduce esta versión

Los perfiles son los exportados y utilizados en la revisión final: Rahue conserva la geometría ligada al outlet del SHP y Curanilahue utiliza su cuenca completa de aproximadamente 119,19 km². Los límites inferiores de ajuste y los comienzos de tiempo se registran por separado. No se cambian las áreas ni se redelinean cuencas al ejecutar la síntesis.

La calibración uniforme utiliza `K_G = (χ_PT − χ_B)/(A₀^θ · T_años)`. Se mantienen los K_P guardados, derivados del estimador de pendiente/tasa costera. La media central utiliza los K de los Pₜ físicos distintos de cada caso. El cierre individual reproduce por construcción la edad asignada a la terraza; no es una validación cronológica independiente.

El modelo compuesto de Reloca mantiene una terraza local MIS 5e. El K del dominio superior procede del promedio de Curanilahue por MIS. El tiempo de ese dominio se resta de los 125 ka locales para calibrar un K inferior positivo en cada alternativa Pₜ. El promedio de esos K conserva residuos que se expresan como RMSE geométrico. C2 es la reconstrucción geológica preferida; C3/C4 tienen menor dispersión y C1 es inadmisible.

Para K uniforme, `τ_c = (χ_KP − χ_c)/(K · A₀^θ · 1000)`, en ka. Para Reloca se integra por dominios y se resta la misma primitiva evaluada en c. Esto conserva la distinción entre el intervalo B–Pₜ que calibra K y el tramo c–KP cuyo tiempo se informa.

U_KP utiliza desnivel y tiempo entre los mismos extremos. Su interpretación como cota superior condicionada supone una pendiente previa mucho menor que la del tramo actual. El término de relieve heredado y la evolución del nivel de base restringen la inferencia tectónica. La incisión por ventanas usa `E = 1000 K_G ksn`; esta versión conserva sus ajustes de entrada y repite la conversión, sin afirmar que volvió a estimar todas las pendientes.

## Niveles de comprobación

1. **Cálculo desde derivados:** K centrales, cierre B–Pₜ individual, 27 escenarios temporales, C2 y otras fronteras, U_KP y conversión de 216 ajustes de incisión. Se cotejan con los resultados finales guardados.
2. **Notebooks de esta carpeta:** se verifica el formato y la ejecución de sus celdas. El notebook 03 abre la distribución pública del modelo 3D; regenerarlo requiere las dependencias y el DEM.
3. **Herramientas interactivas:** se conservan por separado su modelo computacional, interfaz, datos y generador. La selección manual no se sustituye por una regresión automática.
4. **Extracción integral desde DEM:** el archivo histórico contiene los notebooks y módulos utilizados. Esta publicación no vuelve a ejecutar el drenaje D8, los cruces espaciales de terrazas ni los datos privados de geocronología.

## Procedencia de los instrumentos

Los módulos `swath_b_3d.py`, `visor_rios_swath.py` y los tres controladores JavaScript proceden del visor del autor ya incorporado al repositorio en octubre. La adaptación elimina el descubrimiento de rutas antiguas y permite suministrar las entradas por configuración. El HTML público se regenera con los ocho B adoptados y sin el catálogo espacial privado.

El explorador KP conserva el modelo y la interfaz HTML de la revisión del 6 de noviembre. La distribución pública reduce las entradas a 18 perfiles adoptados, sus tributarios y 16 K centrales de referencia. Las variantes P5, K por Pₜ y calibraciones alternativas quedan en el archivo de investigación o en las tablas de sensibilidad, según corresponda. El borrador público se inicia vacío.

La implementación original de extracción topográfica utilizó TopoAnalysis. Su versión o commit exacto debe documentarse desde el entorno de investigación antes de afirmar una reproducción integral. La bibliografía teórica y las referencias de software utilizadas se conservan en la memoria; este paquete no sustituye esas atribuciones.
