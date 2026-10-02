# Geomorfología fluvial y terrazas costeras

Notebooks de Python para analizar topografía, perfiles longitudinales de ríos, χ, pendiente normalizada (ksn), candidatos a knickpoints, terrazas y relaciones espaciales con fallas y unidades geológicas.

Proyecto de investigación en desarrollo de Raimundo Sanchez Mosqueira, estudiante de Geología de la Universidad de Chile. El flujo se desarrolla en el contexto de una memoria de título sobre geomorfología fluvial y evolución del relieve costero de Chile central.

## Estado del proyecto

Los notebooks son herramientas de investigación que requieren configurar rutas, datos y decisiones de análisis para cada aplicación. Todavía no se ha verificado su ejecución completa en un entorno nuevo. Las figuras y resultados guardados en un notebook no reemplazan esa verificación.

Esta documentación describe las versiones inspeccionadas el 2 de octubre de 2026. Los seis notebooks se incluyen en `notebooks/`, sin salidas guardadas. Las rutas personales fueron sustituidas por rutas de ejemplo bajo `C:\datos`; deben adaptarse antes de ejecutar. Los originales de investigación se conservaron sin cambios.

## Flujo de trabajo

| Orden | Notebook | Función |
| --- | --- | --- |
| 1 | `AnalisisTopografico.ipynb` | Carga del DEM, relleno de depresiones, dirección de flujo D8, área de drenaje y delimitación de cuencas desde outlets. |
| 2 | `Perfiles_chi_ksn.ipynb` | Extracción de cauces principales, perfiles longitudinales, integración de χ, ksn local, sensibilidad a ventanas y exportación de perfiles y redes. |
| 3 | `AnalisisPerfiles.ipynb` | Delimitación de tramos, sensibilidad a la concavidad, ajustes χ–elevación y comparación y selección manual de candidatos a knickpoints. |
| 4 | `AnalisisTerrazas.ipynb` | Asociación de terrazas con perfiles, cruces por cota, selección de tramos y escenarios condicionales de K y tiempo de propagación. |
| 5 | `Interpretaciones_K_deltaT.ipynb` | En la versión revisada: carga y revisión de terrazas y knickpoints, lectura de selecciones desde Excel y visualización. El nombre no implica que implemente por sí solo el cálculo completo de K y Δt. |
| Complementario | `Perfiles_y_fallas.ipynb` | Transectas topográficas a fallas, cruces río–falla e integración de perfiles con geología, knickpoints y terrazas. Requiere capas adicionales. |

El flujo intercambia archivos; no es necesario compartir la misma sesión de Python entre notebooks. Sí es necesario que las rutas apunten a productos compatibles de una misma corrida.

## Requisitos

Dependencias identificadas en el código: NumPy, pandas, SciPy, Matplotlib, Rasterio, GeoPandas, Shapely, GDAL/OGR/OSR e IPython. Se requiere además Jupyter para trabajar con los notebooks y un motor compatible para leer Excel, por ejemplo openpyxl.

`TopoAnalysis` se utiliza en el procesamiento topográfico y la extracción de perfiles. `topotoolbox` se importa en la comparación de métodos de `AnalisisPerfiles.ipynb`. Antes de distribuir un entorno instalable deben documentarse sus fuentes y versiones exactas, incluyendo posibles modificaciones locales. No se proporciona todavía un entorno de instalación validado.

## Datos y configuración

Los datos de investigación no se incluyen automáticamente con el código. Consulte `data/README.md` para conocer las entradas necesarias.

1. Prepare el entorno y los datos correspondientes al análisis.
2. Ajuste las rutas de entrada y salida de cada notebook: las versiones publicadas contienen rutas de ejemplo de Windows y referencias a carpetas fechadas.
3. Compruebe el sistema de referencia, las unidades, la extensión y la alineación de las grillas. El caso actual utiliza referencias a UTM 18S (EPSG:32718).
4. Configure los outlets, los identificadores de río y los parámetros científicos. En `Perfiles_chi_ksn.ipynb` se observan A0 = 1 000 000 m², θ de referencia = 0,45, umbral de área = 1 000 000 m² y ventana local inicial = 500 m; son parámetros del caso de estudio, no valores universales.
5. Ejecute las celdas de forma secuencial, comprobando las salidas de cada etapa antes de continuar. Algunas etapas requieren selección manual e interacción gráfica.
6. Conserve los manifiestos, parámetros y selecciones que permitan identificar la corrida usada en los resultados posteriores.

## Interpretación y limitaciones

- Los candidatos automáticos a knickpoints requieren revisión geomorfológica.
- Las asociaciones espaciales entre fallas, terrazas y cambios de pendiente no establecen por sí solas causalidad tectónica.
- K y Δt dependen de las hipótesis del modelo y de las selecciones de terrazas y tramos. El tiempo de propagación modelado no equivale automáticamente a una edad tectónica.
- En el escenario implementado en `AnalisisTerrazas.ipynb`, la tasa de levantamiento debe definirse explícitamente y los supuestos de equilibrio y transferencia de K deben activarse de forma deliberada.
- Las versiones revisadas todavía mezclan referencias a distintas fechas de productos; es necesario armonizarlas antes de reproducir el flujo completo.
- `Perfiles_y_fallas.ipynb` incluye celdas de diagnóstico de lectura del DEM. Su ejecución completa no se ha validado en esta revisión.

## Próximos pasos

- [ ] Sustituir rutas personales por una configuración portable.
- [ ] Documentar y fijar el entorno que efectivamente funciona.
- [ ] Incorporar un conjunto pequeño de datos con permiso de redistribución.
- [ ] Verificar el flujo desde un entorno nuevo.
- [ ] Completar las referencias bibliográficas y las versiones de las dependencias científicas.
- [x] Agregar licencia MIT para las contribuciones propias.
- [ ] Preparar una primera versión estable y sus instrucciones de citación.

## Atribución, licencia y contribuciones

El flujo utiliza software y métodos de terceros que deben citarse por separado. TopoAnalysis y TopoToolbox no son desarrollos propios de este proyecto. Las licencias y atribuciones de código o datos externos deben conservarse.

Las contribuciones propias de este proyecto se distribuyen bajo la licencia MIT: consulte `LICENSE`. Las dependencias y los materiales de terceros conservan sus propias licencias. Esta licencia no concede derechos sobre datos externos.

Para informar problemas, abra un issue que incluya el notebook, el mensaje de error, el entorno y los pasos necesarios para reproducirlo. Las propuestas de cambios deben explicar su efecto en el método o en los resultados.

