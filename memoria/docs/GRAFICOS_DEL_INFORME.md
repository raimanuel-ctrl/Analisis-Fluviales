# Reproducción de los gráficos del informe

La revisión toma como referencia las 56 figuras de la memoria final limpia reparada del 7 de noviembre de 2026. El generador [graficos_informe.py](../src/graficos_informe.py) vuelve a dibujar las series y los paneles de gráficos a partir de sus entradas numéricas. Conserva las posiciones adoptadas, los estados potenciales, los estimadores de K y los orígenes temporales. La distribución y el estilo pueden variar respecto del Word; las salidas son gráficos recalculables, no capturas del documento.

La revisión comprende **47 gráficos y esquemas propios**, **una reproducción bibliográfica** y **ocho figuras cartográficas excluidas**. Los esquemas 4–7 utilizan geometría ilustrativa expresamente identificada. En la figura mixta 14 solo se genera el perfil inferior; sus mapas quedan fuera.

## Ejecutar

Desde la carpeta `memoria`:

```sh
python -m src.graficos_informe
```

Cada ejecución crea una carpeta nueva en `resultados/`, con los gráficos numerados en PNG y SVG, `cobertura.json` y las huellas de sus entradas. Para elegir figuras o formatos:

```sh
python -m src.graficos_informe --figuras 11 15 20 34 41 48
python -m src.graficos_informe --figuras 15 20 --formatos png svg pdf
```

El [notebook 09](../notebooks/09_graficos_del_informe.ipynb) permite revisar el inventario y ejecutar el generador. Las entradas adicionales están descritas en [datos/graficos](../datos/graficos/README.md). No se ejecuta cartografía ni se vuelven a delinear cuencas.

## Cobertura por figura

| Figuras | Contenido | Tratamiento |
| --- | --- | --- |
| 1, 3, 12, 17, 25, 27, 32, 40 | Mapas y composiciones cartográficas | Excluidas. Se mantienen como referencia para las composiciones del autor. |
| 2 | Figura 12 de Jara-Muñoz et al. (2015) | Reproducción bibliográfica desde una imagen local autorizada; no se presenta como resultado calculado por este código. |
| 4–7 | Perfiles, propagación, estimadores y procedimiento conceptuales | Redibujo mediante código; escalas ilustrativas, sin datos medidos. |
| 8 | Cuatro ejemplos de cierres sobre perfiles costeros | Perfiles costeros propios, origen del cauce y KP guardados. |
| 9 | Tramo temporal de Curanilahue | Perfil completo, comienzo c y KP adoptado. |
| 10 | Alternativas Pₜ de Pinotalca | Las tres etiquetas originales y sus cruces físicos sobre H. |
| 11 | Los 18 perfiles χ–elevación | Elevación relativa al outlet y traslado gráfico de 650 m entre curvas. |
| 13 | Calibración geométrica en Chanco y Curanilahue | Perfiles, H, B, Pₜ intermedio y Δχ de cada MIS. |
| 14 | Figura mixta de Curanilahue | Solo perfil; mapas a y b excluidos. |
| 15 | K_G y K_P por latitud | Medias por punto físico, separadas por MIS y método. |
| 16 | Comparación de K con referencias externas | Calibraciones propias y serie externa suministrada localmente. |
| 18 | Reloca C2 en elevación–τ | Dos MIS del donante, misma terraza local y origen c. |
| 19 | C2–C4 y residuos de calibración | Tiempos y RMSE geométrico recalculados. |
| 20 | Tiempos por latitud, MIS y método | Los 27 escenarios finales. |
| 21 | Altura de KP por latitud y tiempo | R01 entra en altura/latitud y queda fuera del panel temporal. |
| 22 | Efecto de H en los tiempos del sur | Escenarios DEM/altura publicada guardados. |
| 23 | Perfiles completos y adoptados de R02–R04 | Geometrías del ensayo conservadas por separado. |
| 24 | Tres donantes y KP alternativos del norte | Resultados guardados del ensayo, diferenciados del catálogo final. |
| 26 | Perfil de Rahue y cruces de Pelluhue | Perfil y distancias de intersección guardadas, sin dibujar un mapa. |
| 28 | Incisión en la ventana de 2 km | Ajustes de referencia, con cobertura parcial y KP potenciales identificados. |
| 29 | U_KP frente a tasas de terrazas | Desnivel/tiempo final y serie costera aportada localmente. |
| 30 | Perfil y fronteras de Reloca | Perfil completo y detalle de C1–C4. |
| 31 | Síntesis de los seis receptores del sur | Rango, promedio y mediana por MIS/método. |
| 33 | Sensibilidad a cinco decisiones | Contrastes de KP, recorrido, altura, MIS y método. |
| 34 | Los perfiles elevación–τ en cuatro paneles | Misma integración temporal que los resultados finales; Reloca usa C2. |
| 35 | Diagnóstico de concavidad | Curvas históricas guardadas, con el mismo hash de la fuente utilizada en la figura final. |
| 36–38 | Perfiles de R01, R07 y R17 | R07 permanece sin KP y R17 conserva su condición potencial. |
| 39 | OD y alícuotas de doce muestras | Tabla 1 publicada de Jara-Muñoz et al. (2015); sin usar Ages_Final. |
| 41–43 | Sensibilidad de R03, Reloca y Rahue | Perfiles, rectas, incisión y cobertura de ocho ventanas. |
| 44 | Sensibilidad de incisión a la ventana | 216 combinaciones; color saturado en ±100 %, datos completos en CSV. |
| 45 | Reloca uniforme frente a C2 | Misma posición de KP y mismo inicio temporal. |
| 46–47 | Buchupureo y Curanilahue en elevación–τ | K_G/K_P y MIS separados, conservando el tramo anterior a c como contexto. |
| 48–56 | Regresiones completas y de 3 km | Coeficientes guardados, exclusión inferior, referencia manual y KP adoptado. Los cruces exploratorios se muestran solo dentro del dominio admitido. |

El inventario completo, con los pies de figura del documento, está en [inventario_figuras.json](../datos/graficos/inventario_figuras.json). La verificación de la ejecución está en [VERIFICACION_GRAFICOS.json](../VERIFICACION_GRAFICOS.json).

## Comparaciones que requieren entradas locales

Tres figuras conservan una dependencia externa explícita:

| Figura | Archivo local | Campos o contenido |
| --- | --- | --- |
| 2 | `Jara_Munoz_2015_Fig12.png` | Imagen de la figura bibliográfica con su atribución. |
| 16 | `comparacion_K_externa.csv` | `sector`, `K`, `error_K`; valores de la comparación externa. |
| 29 | `terrazas_comparacion.csv` | `latitud_sur`, `U_mm_anio`, `error_U_mm_anio`; serie costera utilizada en la comparación. |

Los archivos se colocan en `entradas_privadas/`, excluida de Git, o se indica otra carpeta:

```sh
python -m src.graficos_informe --entradas-privadas RUTA_LOCAL
```

La copia local de trabajo contiene las entradas necesarias y permitió comprobar las tres figuras. La distribución pública mantiene sus contratos de entrada; si faltan, `cobertura.json` identifica exactamente las figuras pendientes y el resto de los gráficos se genera normalmente. No se sustituyen las series externas por valores de ejemplo.

## Qué se verificó

Se ejecutó el conjunto completo sin generar los ocho mapas. Además de comprobar las imágenes y los SVG, se reconstruyeron por mínimos cuadrados los **68 ajustes** de las regresiones completa/local y los **216 ajustes de incisión** sobre sus nodos originales. Coincidieron los números de nodos y los coeficientes guardados, dentro de la precisión numérica. La fuente de concavidad conserva el SHA-256 registrado para la figura 35 final.

Esta comprobación verifica la correspondencia entre código, entradas y cantidades representadas. La aceptación geomorfológica de B, Pₜ y KP y la interpretación de K y U_KP conservan los supuestos expuestos en el README y en la memoria.
