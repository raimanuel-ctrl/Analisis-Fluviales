# Tablas de presentación de la memoria

Transcripción de las tablas seleccionadas de la memoria final limpia reparada del 7 de noviembre de 2026. Los nombres de archivo describen su contenido y conservan la numeración del documento. Las tablas 15 y 19 añaden la columna `ID KP` con las etiquetas públicas, sin cambiar los datos originales. `indice_tablas.csv` permite relacionarlas sin depender del nombre de una carpeta de trabajo.

| Tabla | Archivo | Contenido |
| --- | --- | --- |
| 8 | [Tabla_08_alturas_dem_y_publicadas.csv](Tabla_08_alturas_dem_y_publicadas.csv) | Alturas y diferencias DEM/catálogo. |
| 9 | [Tabla_09_k_por_interseccion.csv](Tabla_09_k_por_interseccion.csv) | K_G de cada etiqueta Pₜ y rango nominal. |
| 10 | [Tabla_10_k_centrales.csv](Tabla_10_k_centrales.csv) | K_G y K_P centrales por río/MIS. |
| 11 | [Tabla_11_sensibilidad_altura_k.csv](Tabla_11_sensibilidad_altura_k.csv) | Sensibilidad de K_G a H. |
| 13 | [Tabla_13_k_y_tiempos_calibracion_propia.csv](Tabla_13_k_y_tiempos_calibracion_propia.csv) | K y τ en ríos con calibración propia. |
| 14 | [Tabla_14_reloca_fronteras.csv](Tabla_14_reloca_fronteras.csv) | C1–C4, tiempos y residuos. |
| 15 | [Tabla_15_knickpoints_y_tiempos.csv](Tabla_15_knickpoints_y_tiempos.csv) | Los 18 cauces, altura de KP, tiempos y estado. |
| 16 | [Tabla_16_sensibilidad_altura_tiempos.csv](Tabla_16_sensibilidad_altura_tiempos.csv) | Efecto de H en los tiempos. |
| 17 | [Tabla_17_reloca_sensibilidad_alturas.csv](Tabla_17_reloca_sensibilidad_alturas.csv) | Altura y donante en Reloca C2. |
| 18 | [Tabla_18_incision_y_tasas_terrazas.csv](Tabla_18_incision_y_tasas_terrazas.csv) | Incisión B–Pₜ y tasas costeras. |
| 19 | [Tabla_19_alzamiento_medio_tramo.csv](Tabla_19_alzamiento_medio_tramo.csv) | Desnivel, τ_c y U_KP. |
| 20 | [Tabla_20_resumen_seis_receptores.csv](Tabla_20_resumen_seis_receptores.csv) | Estadísticos del sur por método/MIS. |
| 33 | [Tabla_33_origen_temporal_y_recorrido.csv](Tabla_33_origen_temporal_y_recorrido.csv) | Referencias temporales y sensibilidad al recorrido. |

Estas tablas conservan la coma decimal y el redondeo del documento. Los valores ausentes aparecen como raya, igual que en la memoria. Para calcular se utilizan exclusivamente los CSV numéricos de [datos](../datos/README.md), donde un dato numérico ausente se representa por un campo vacío.

Las unidades de K de las tablas 9–14 son ×10⁻⁵ m⁰·²/año; las edades y tiempos están en ka y las tasas en mm/año. Las alturas y distancias verticales se expresan en m; la tabla 33 identifica expresamente s_c en km. La condición potencial se conserva mediante su columna o asterisco original. `Tabla_15` enumera los 18 ríos, incluidos R01 sin K y R07 sin KP.
