# Análisis fluvial y síntesis de resultados

Versión pública del 22 de octubre de 2026. Siete notebooks con módulos de apoyo, sin salidas guardadas, rutas personales ni datos de colaboración no publicados. Los originales de trabajo permanecen sin modificar.

## Empezar por la síntesis reproducible

Abrir [07_Sintesis_resultados_20261022.ipynb](Notebooks_principales_20261019/07_Sintesis_resultados_20261022.ipynb). Incluye sus CSV derivados en `datos_sintesis`, con huellas SHA-256. No necesita una aplicación interactiva, documentos privados ni datos externos de comparación. También puede ejecutarse desde esta carpeta:

```text
python -m pip install -r requirements_sintesis.txt
python sintesis_20261022.py
```

Cada ejecución crea una carpeta nueva en `Productos/sintesis_20261022/`, con CSV, PNG, SVG y verificación. Los productos locales no se versionan. La prueba reproduce las tablas de tiempos y tasas usadas en la memoria del 22 de octubre. **No vuelve a extraer la red desde el DEM ni valida de manera independiente la interpretación tectónica.**

## Orden de trabajo

| Orden | Notebook | Función y entradas |
| --- | --- | --- |
| 01 | [Análisis topográfico](Notebooks_principales_20261019/01_Analisis_Topografico.ipynb) | DEM, D8, áreas y perfiles originales; revisar configuración y TopoAnalysis externo. |
| 02 | [Perfiles](Notebooks_principales_20261019/02_Perfiles.ipynb) | Morfometría, KP, litología y terrazas descriptivas; necesita productos y manifiestos de 01 y entradas cartográficas. |
| 03 | [K U y τ](Notebooks_principales_20261019/03_K_U_tau_20261019.ipynb) | Calibraciones B–Pₜ, alternativas y cronología MIS; necesita decisiones manuales y los productos declarados. Las comparaciones externas están desactivadas y sus valores no se distribuyen. |
| 04 | [ksn y fallas](Notebooks_principales_20261019/04_ksn_Fallas.ipynb) | Cauces seleccionados y comparación estructural; no recalibra K. |
| 05 | [Red completa](Notebooks_principales_20261019/05_Red_completa_fallas_20261019.ipynb) | Red del polígono y cauces que cruzan fallas; necesita D8 y área del DEM completo, polígono y trazas. |
| 06 | [Figuras](Notebooks_principales_20261019/Figuras.ipynb) | Generador general y cuatro fronteras de Reloca; necesita sus productos y cartografía externos. No sustituye la síntesis C2 de 07. |
| 07 | [Síntesis actual](Notebooks_principales_20261019/07_Sintesis_resultados_20261022.ipynb) | K promedio por punto físico, τ con C2, comparación método/MIS, residuos de calibración, tasas bajo KP y figuras estáticas. Incluye sus entradas derivadas. |

01–06 son copias del conjunto principal vigente localizado, con adaptaciones de publicación; **no se reejecutaron**. Se comprueban esquema y sintaxis por separado de la prueba numérica de 07. Los visores de revisión manual son opcionales y no constituyen el mecanismo de reproducción de la síntesis.

Iniciar Jupyter desde `actual_20261022` o `Notebooks_principales_20261019`. Reiniciar el kernel entre notebooks: algunos módulos tienen nombres compartidos pero versiones específicas por rama. `python publicacion_config.py` lista entradas disponibles y externas. No se descargan datos ni instalan paquetes al abrir un notebook. La reproducción integral desde el DEM sigue requiriendo insumos externos y un entorno validado para esas etapas.

## Convenciones de la versión actual

- θ = m = 0,40, n = 1 y A₀ = 10⁶ m². χ se integra trapezoidalmente; `ksn = Mχ A₀^θ`; K en m⁰·²/año.
- K_G representa propagación y K_P, pendiente y tasa de terraza. Se conservan identificadores heredados para no romper la interfaz original.
- Ocho casos río/MIS y 24 etiquetas Pₜ corresponden a 18 puntos físicos distintos. Cada punto físico contribuye una vez al promedio; primero se promedia K y después se calcula τ. No se mezclan MIS.
- Reloca utiliza C2: K local 5e abajo y Curanilahue 5c o 5e arriba. Los residuos muestran que no se conservan simultáneamente las tres referencias inferiores. C3 y C4 quedan como controles.
- R02–R03 reciben R04–5e; R10–R18 reciben R09 por MIS. No son dataciones independientes. R01 no tiene K asignado y R07 no tiene un KP aceptado.
- R09 conserva el escenario de área modificada y K sin recalibrar; Rahue comienza en B N07. `Geometria_origenes.csv` declara cada origen. Dibujar el perfil hasta la costa no cambia el origen temporal.
- Las tasas usan hasta 2.000 m inmediatamente bajo el KP, no B–Pₜ. Se reportan longitud y R². `1.000 K_G ksn` es incisión modelada en mm/año; sólo equivale a U bajo equilibrio local.
- El gráfico de tasas no superpone valores costeros publicados. El 3D mira hacia N10°E sin hillshade.
- Los rangos de escenarios no son intervalos de confianza; las correlaciones MIS no son nuevas dataciones del cauce. No se confirma un pulso único de 0,8–1 Ma ni la inactividad de las fallas.

## Datos y permisos

`datos_sintesis` contiene geometría derivada, KP proyectados, calibraciones propias, orígenes y tablas de referencia seleccionadas. La cronología conserva las correlaciones de Jara-Muñoz et al. (2015), DOI [10.1016/j.quascirev.2015.01.005](https://doi.org/10.1016/j.quascirev.2015.01.005), no sus mediciones primarias de luminiscencia. CRS EPSG:32718; el datum vertical absoluto del DEM requiere comprobación.

El DEM y la compilación geológica de la publicación anterior están en [data](../../data/), accesibles desde las rutas configuradas; el DEM utiliza Git LFS. Los datos fuente no heredan automáticamente la licencia del código. No se publican correos, manuscritos, libros privados, notebooks de colaboradores ni sus coeficientes no publicados. Las comparaciones que los requieren no se reemplazan por valores inventados.

Para ejecutar 01–06 se necesitan también fallas, polígonos, catálogos de terrazas, decisiones B/KP y manifiestos compatibles. Las fuentes derivadas de 07 no sustituyen esos insumos. La posibilidad de repetir el cálculo tampoco valida automáticamente la selección geomorfológica manual.

## Procedencia y comprobaciones

`procedencia_20261022.json` registra originales y transformaciones; `manifest_publicacion_20261022.json`, los archivos finales. Ver [REVIEW_20261022.md](../../docs/REVIEW_20261022.md). `requirements_sintesis.txt` cubre la síntesis mínima; la lista de investigación no promete un entorno completo validado.
