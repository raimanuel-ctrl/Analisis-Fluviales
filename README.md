# Perfiles fluviales y evolución del relieve costero de Chile central

Material computacional de la memoria de título en Geología de Raimundo Sánchez Mosqueira, Universidad de Chile. Incluye los notebooks, instrumentos de selección y datos derivados utilizados para comparar los perfiles de 18 cauces entre 35°S y 36,5°S.

**Versión de entrega: 2026.11.09.** Los resultados corresponden a la memoria final del 7 de noviembre. La carpeta [memoria](memoria/README.md) es la única entrada al flujo vigente.

## Accesos principales

| Material | Acceso |
| --- | --- |
| Ocho notebooks ordenados | [Ruta de trabajo 01–08](memoria/README.md#empezar) |
| Reproducción de resultados | [08 · Síntesis y tablas](memoria/notebooks/08_sintesis_y_tablas_memoria.ipynb) |
| Catálogo de 17 KP | [Posiciones adoptadas y potenciales](memoria/datos/catalogo_kp.csv) |
| Correspondencia de identificadores | [Códigos públicos y originales](memoria/datos/equivalencias_identificadores_kp.csv) |
| Selección geomorfológica de B | [Notebook 03 y visor 3D](memoria/notebooks/03_seleccion_B_y_visor_3d.ipynb) |
| Explorador interactivo de perfiles/KP | [HTML autónomo](memoria/visores/Explorador_KP.html) |
| Relieve 3D de los B adoptados | [HTML autónomo](memoria/visores/Visor_B_3D.html) |
| Tablas de la memoria | [Índice y trece tablas finales](memoria/tablas_memoria/README.md) |
| Datos y unidades | [Diccionario](memoria/docs/DICCIONARIO_DATOS.md) |
| Qué se reprodujo | [Metodología y comprobaciones](memoria/docs/METODOLOGIA_Y_REPRODUCIBILIDAD.md) |
| Referencia para el manuscrito | [Disponibilidad y cita del código](memoria/docs/REPOSITORIO_EN_LA_MEMORIA.md) |

Para abrir una aplicación, descargar el HTML y abrirlo en un navegador. GitHub muestra el código del archivo; no ejecuta el visor en su vista de repositorio. La distribución funciona localmente sin enviar datos a un servidor.

## Reproducir la síntesis

```sh
git clone https://github.com/raimanuel-ctrl/Analisis-Fluviales.git
cd Analisis-Fluviales/memoria
python -m pip install -r requirements.txt
python -m src.analisis
```

Cada ejecución crea una carpeta nueva con seis tablas numéricas, tres figuras de control en PNG/SVG y la verificación de resultados. La síntesis utiliza datos derivados incluidos. El DEM es necesario para regenerar el relieve 3D y para las etapas topográficas originales; permanece en `data/dem/` mediante Git LFS.

## Alcance científico

Se utilizan θ = 0,40, n = 1 y A₀ = 10⁶ m². Las ocho asociaciones río/MIS reúnen 18 Pₜ físicos distintos. El promedio central de K cuenta cada punto una vez y conserva separados MIS 5c y 5e. Reloca utiliza C2 con dos dominios efectivos; los receptores del sur reciben el K de la cuenca completa de Curanilahue, con MIS 5e como escenario principal.

Los tiempos de propagación se calculan desde los comienzos de tramo adoptados. U_KP = (z_KP − z_c)/τ_c se interpreta como una cota superior condicionada del alzamiento medio bajo una pendiente inicial mucho menor que la actual. Las tasas de incisión, las tasas costeras y U_KP conservan sus definiciones y duraciones de integración. La reproducción del cálculo no demuestra por sí sola una edad tectónica, un pulso regional único ni la inactividad de una falla.

R01 conserva KP sin K adoptado; R07 tiene K y carece de KP aceptado. R11, R14 y R17 permanecen potenciales. Los seis receptores principales del sur son R10, R12, R13, R15, R16 y R18.

## Organización y conservación

- [memoria/](memoria/README.md): versión vigente, datos derivados, herramientas y tablas finales.
- [data/](data/README.md): DEM, drenaje y geología publicados previamente, con sus condiciones de origen.
- [archivo/](archivo/README.md): versiones de octubre y antecedentes metodológicos, conservados para trazabilidad.

Las versiones anteriores quedan identificadas como históricas para evitar mezclar K, geometrías y orígenes temporales. Los originales de trabajo se conservan fuera de esta reorganización. Las [notas de selección](memoria/docs/SELECCION_DE_CONTENIDOS.md) explican qué material se mantuvo activo.

El código propio está cubierto por [MIT](LICENSE). Las fuentes científicas conservan su autoría y condiciones. Los catálogos espaciales privados de terrazas, Ages_Final, correspondencia, manuscritos y valores no publicados de colaboradores quedan fuera de la distribución. La reproducción integral desde el DEM requiere insumos adicionales y el entorno de las etapas originales.

## Citar

Sánchez Mosqueira, R. (2026). *Analisis-Fluviales: material computacional de la memoria de geomorfología fluvial* (versión 2026.11.09) [Código y datos derivados]. GitHub. https://github.com/raimanuel-ctrl/Analisis-Fluviales

Los metadatos están en [CITATION.cff](CITATION.cff). Para la entrega de la memoria se debe conservar también el identificador del commit utilizado.
