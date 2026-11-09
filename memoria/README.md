# Material computacional de la memoria

Versión preparada el 9 de noviembre de 2026 a partir de la memoria final del 7 de noviembre. Esta carpeta reúne los notebooks, las selecciones manuales y los datos derivados que sostienen sus resultados. El material de investigación anterior se conserva en [archivo](../archivo/README.md).

## Empezar

Para reproducir los resultados finales, ejecutar [08 · Síntesis y tablas](notebooks/08_sintesis_y_tablas_memoria.ipynb). Para revisar el razonamiento completo, seguir el orden 01–08. Cada notebook tiene una función y utiliza las mismas tablas de entrada.

```sh
cd memoria
python -m pip install -r requirements.txt
python -m src.analisis
```

La ejecución crea una carpeta nueva dentro de `resultados/`. Recalcula los K centrales por Pₜ físico, los tiempos, el modelo compuesto de Reloca, el alzamiento del tramo y la conversión de los ajustes de incisión guardados. Genera seis tablas, tres figuras de síntesis en PNG/SVG y un registro de verificación. Los archivos generados permanecen locales.

| Orden | Notebook | Función |
| --- | --- | --- |
| 01 | [Perfiles y control topográfico](notebooks/01_perfiles_y_control_topografico.ipynb) | Revisar los 18 perfiles derivados y sus coordenadas. |
| 02 | [Catálogo e identificación de KP](notebooks/02_catalogo_e_identificacion_kp.ipynb) | Revisar las posiciones adoptadas, los potenciales y sus identificadores. |
| 03 | [Selección B y visor 3D](notebooks/03_seleccion_B_y_visor_3d.ipynb) | Documentar la revisión geomorfológica de B y abrir o regenerar el relieve. |
| 04 | [Calibración de K](notebooks/04_calibracion_k_y_referencias_mis.ipynb) | Calcular medias por Pₜ físico y comprobar el cierre individual B–Pₜ. |
| 05 | [Reloca C2](notebooks/05_reloca_modelo_compuesto_C2.ipynb) | Recalcular dominios inferior/superior, alternativas C1–C4 y residuos. |
| 06 | [Tiempos de propagación](notebooks/06_tiempos_de_propagacion.ipynb) | Integrar tiempos desde los orígenes temporales adoptados. |
| 07 | [Alzamiento e incisión](notebooks/07_alzamiento_e_incision.ipynb) | Calcular U_KP, conservar el contraste de cota absoluta y evaluar los ajustes por ventana. |
| 08 | [Síntesis y tablas](notebooks/08_sintesis_y_tablas_memoria.ipynb) | Ejecutar la reproducción y cotejar las tablas de la memoria. |

## Herramientas interactivas

- [Explorador de KP](visores/Explorador_KP.html): archivo HTML autónomo con perfiles, mapa, tributarios de área ≥1 km² y regresiones exploratorias. La distribución abre con un borrador vacío. Los puntos adoptados se identifican mediante `KP-R##`; los clics nuevos permanecen exploratorios.
- [Visor de B en 3D](visores/Visor_B_3D.html): relieve y bandas N–S/O–E para los ocho B adoptados. Conserva las funciones de giro, dimensiones, exageración vertical y selección manual del visor utilizado en la investigación. La versión pública omite las ubicaciones originales del catálogo de terrazas.

Descargar los HTML para utilizarlos en el navegador; la vista de archivos de GitHub no ejecuta aplicaciones. Ambos funcionan sin servidor. Para regenerarlos desde esta carpeta:

```sh
python -m pip install -r requirements_visores.txt
python -m src.visores kp
python -m src.visores b3d
```

El visor 3D requiere el DEM de `data/dem/`, descargado mediante Git LFS. Su ruta se configura en [configuracion.json](configuracion.json). Las referencias de terrazas pueden añadirse desde un CSV autorizado localmente. El HTML generado con ese catálogo incorpora sus puntos y debe conservarse localmente.

El explorador calcula un K por medición. El modelo compuesto C2 se reproduce en los notebooks y no se aplica automáticamente en esa aplicación.

## Datos y tablas

[datos](datos/README.md) contiene los perfiles, el catálogo de KP, las selecciones B/Pₜ, las calibraciones y los escenarios con precisión numérica. [tablas_memoria](tablas_memoria/README.md) reúne trece tablas transcritas de la versión final, con su numeración y redondeo de presentación. El [diccionario de datos](docs/DICCIONARIO_DATOS.md) define los campos y las unidades. Las tablas redondeadas sirven para cotejar el texto; los cálculos leen los datos de precisión completa.

La reproducción parte de datos derivados. La extracción D8 original desde el DEM y las asociaciones espaciales con terrazas requieren sus insumos y su entorno original, conservados o identificados en el archivo metodológico. La comprobación numérica de esta carpeta no constituye una nueva validación geomorfológica de las selecciones.

## Convenciones científicas

θ = 0,40; n = 1; A₀ = 10⁶ m²; K en m⁰·²/año. K_G designa propagación y K_P, pendiente y tasa costera. Se promedia K sobre Pₜ físicos distintos dentro de cada caso río/MIS, antes de calcular τ. El caso principal del sur transfiere K_G de Curanilahue–5e, con su cuenca completa. MIS 5c permanece como contraste.

Los tiempos parten de c, el comienzo del tramo temporal adoptado; B delimita la calibración con Pₜ. Estas referencias se guardan por separado. U_KP = (z_KP − z_c)/τ_c se interpreta como una cota superior condicionada del alzamiento medio bajo una pendiente inicial mucho menor que la actual. K_G ksn representa incisión modelada. Las edades y las tasas conservan los supuestos del modelo y el datum vertical absoluto del DEM sigue sin confirmación.

R01 conserva un KP sin K adoptado. R07 tiene calibraciones de K y carece de KP aceptado. R11, R14 y R17 son potenciales. La síntesis principal del sur incluye R10, R12, R13, R15, R16 y R18.

Ver [reproducibilidad](docs/METODOLOGIA_Y_REPRODUCIBILIDAD.md), [selección de contenidos](docs/SELECCION_DE_CONTENIDOS.md) y [cómo referenciar el repositorio en la memoria](docs/REPOSITORIO_EN_LA_MEMORIA.md).
