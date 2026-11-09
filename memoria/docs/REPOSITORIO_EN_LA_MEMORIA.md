# Referencia y descripción para la memoria

Texto propuesto para metodología o disponibilidad del código:

> Los notebooks y módulos de procesamiento se organizaron en el repositorio Analisis-Fluviales (Sánchez Mosqueira, 2026). La versión de la memoria reúne los perfiles derivados de los 18 cauces, las selecciones de B y Pₜ, el catálogo de knickpoints y los parámetros necesarios para reproducir los coeficientes centrales de K, los tiempos de propagación y las estimaciones de alzamiento medio del tramo. Los identificadores de presentación conservan una correspondencia con los códigos originales de las selecciones.
>
> La distribución incluye el código del visor 3D empleado en la revisión de B y del explorador de perfiles y KP. Las tablas finales se identifican por su número y contenido. La reproducción numérica parte de datos derivados y se distingue del procesamiento integral del DEM, que requiere el entorno y los insumos de sus etapas originales. Los catálogos privados y los datos no publicados de colaboradores se mantienen fuera de la distribución.

Referencia bibliográfica propuesta:

Sánchez Mosqueira, R. (2026). *Analisis-Fluviales: material computacional de la memoria de geomorfología fluvial* (versión 2026.11.09) [Código y datos derivados]. GitHub. https://github.com/raimanuel-ctrl/Analisis-Fluviales

Para fijar la versión de entrega, utilizar además el identificador del commit que contiene esta actualización. `CITATION.cff` proporciona los metadatos bibliográficos del código. No se asigna un DOI inexistente.

## Nombres de tablas e identificadores

Las tablas públicas conservan su numeración en el Word final del 7 de noviembre. Su correspondencia está en `tablas_memoria/indice_tablas.csv`. Si cambia la numeración en la edición definitiva, debe actualizarse ese índice y el nombre de presentación, manteniendo el dato fuente.

En las tablas de KP se puede presentar `KP-R01`, `KP-R02`, etc., acompañados por el río y una columna de estado. La condición potencial de R11/R14/R17 se expresa en esa columna, independientemente del identificador. R07 figura como cauce sin KP aceptado; no se le asigna un punto por completar la secuencia. Las figuras o anexos antiguos se pueden relacionar mediante `equivalencias_identificadores_kp.csv`.

Esta propuesta se entrega separada del manuscrito para que la referencia pueda insertarse en la ubicación definitiva, sin modificar el archivo original de la memoria.
