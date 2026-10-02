# Datos requeridos

Los datos científicos no están incluidos en esta distribución.

| Etapa | Entradas |
| --- | --- |
| Topografía | DEM georreferenciado y coordenadas de outlets; el caso actual referencia `TDM_cut_UTM18.tif`. |
| Perfiles | DEM original y productos de área y dirección de flujo. |
| Análisis de perfiles | CSV de perfiles y `manifest.json`. |
| Terrazas | Perfiles, `catalogo_candidatos.csv`, `Terrazas_revisadas_UTM18S.xlsx` y `Revision_muestras_Jara_2015.xlsx`. |
| Selecciones | `terrazas_kp_por_rio.xlsx` y entradas anteriores. |
| Fallas y geología | DEM, capas de fallas, red fluvial, unidades geológicas, terrazas y knickpoints. |

Revisar las columnas y validaciones en cada notebook. Los shapefiles requieren sus archivos asociados.

Para cada conjunto de datos, documentar fuente, autoría, licencia, CRS, unidades y transformaciones. Distribuir únicamente datos con permiso; para los demás, indicar cómo obtenerlos.

Estructura propuesta: `raw/` para datos originales externos y `processed/` para productos grandes, ambas excluidas por `.gitignore`; `example/` para un futuro ejemplo pequeño redistribuible.

Las selecciones manuales y parámetros deben conservarse en archivos pequeños versionados, por ejemplo en `config/`. No tratarlos como productos descartables.
