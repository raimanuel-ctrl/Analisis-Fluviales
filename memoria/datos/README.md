# Datos derivados de la versión final

| Archivo | Unidad de registro | Uso |
| --- | --- | --- |
| `catalogo_rios.csv` | 18 cauces | Nombres, estado de KP, disponibilidad de K y dimensiones del perfil exportado. |
| `perfiles/R##_perfil.csv` | Nodos de un cauce | Coordenadas, elevación, distancia, área y χ guardados. |
| `catalogo_kp.csv` | 17 KP adoptados o potenciales | Posición y método de localización; intervalos originales de los ajustes. |
| `equivalencias_identificadores_kp.csv` | Alias de KP | Correspondencia entre etiqueta pública y códigos actuales o previos. |
| `selecciones_b.csv` | 8 asociaciones B/MIS | Centros del visor y geometría de calibración. |
| `calibraciones_k_por_pt.csv` | 24 etiquetas Pₜ | 18 puntos físicos distintos, K_G/K_P, geometría y alturas. |
| `referencias_mis.csv` | 8 referencias genéricas | Edad MIS y fuente publicada, sin ubicaciones de terrazas. |
| `origenes_temporales.csv` | 18 cauces | Inicio del tiempo y corte del ajuste como referencias separadas. |
| `escenarios_tiempo.csv` | 27 escenarios río/MIS | K, posición del KP, origen del tiempo y valores finales de referencia. |
| `reloca_fronteras.csv` | 4 fronteras | Geometría C1–C4. |
| `reloca_sensibilidad_referencia.csv` | 8 frontera/MIS | Tiempos y residuos guardados para cotejar la recalibración. |
| `ajustes_incision_ventanas.csv` | 216 ajustes guardados | Ventanas, ksn, K, cobertura y tasa de referencia. |
| `sensibilidad_altura_calibracion.csv` | Alternativas de H/Pₜ | Contraste DEM/altura publicada, conservando sus parámetros. |
| `procedencia.json` | Versión | Parámetros, alcance y huellas de las fuentes. |

Los CSV usan coma como separador, punto decimal y UTF-8. Las coordenadas propias del cauce, KP, B y Pₜ son productos derivados en EPSG:32718. Un Pₜ es una intersección calculada sobre el cauce. Estas posiciones no son las ubicaciones medidas del catálogo original de terrazas, que queda excluido.

Las referencias cronológicas corresponden a 106 ±5 ka para MIS 5c y 125 ±5 ka para MIS 5e, según Jara-Muñoz et al. (2015), tabla 2, p. 179, DOI [10.1016/j.quascirev.2015.01.005](https://doi.org/10.1016/j.quascirev.2015.01.005). La asociación río/nivel es la correlación adoptada en esta memoria. Las edades de los niveles son entradas de calibración.

Las etiquetas `primero`, `intermedio_real` y `ultimo` pueden designar el mismo Pₜ físico. `id_pt` permite contarlo una vez. Las posiciones previas del KP se conservan como alias de trazabilidad, no como selecciones adicionales vigentes.

El código propio se distribuye con MIT. La presencia de un archivo derivado no modifica los derechos ni las condiciones de sus fuentes científicas. No se distribuyen Ages_Final, correspondencia, manuscritos ni coeficientes externos no publicados.
