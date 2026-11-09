# Convenciones e identificadores

`rid` es el número del río, de 1 a 18; `id_rio` lo presenta como `R01`–`R18`. Los nombres del río se guardan en `catalogo_rios.csv`.

`id_kp` identifica la selección actual como `KP-R##`. Es una etiqueta estable de presentación; no contiene fecha, método ni estado. `estado_kp` distingue `reconocido`, `potencial` y, en el catálogo de ríos, `sin_kp_aceptado`. Así se puede cambiar una interpretación sin alterar artificialmente la identidad del rasgo. Cada cambio de posición debe registrarse como una revisión con su correspondencia. `id_fuente_actual` mantiene el código anterior de la selección adoptada.

`equivalencias_identificadores_kp.csv` relaciona la etiqueta pública con los códigos históricos. `rol=seleccion_adoptada` corresponde al identificador de la posición actual; `referencia_previa` identifica antecedentes de la revisión. La tabla permite rastrear nombres, no equiparar coordenadas de selecciones diferentes.

`id_b`, por ejemplo `B-R09-5e`, identifica una asociación de calibración del cauce con un MIS. `id_pt`, por ejemplo `PT-R09-5e-01`, identifica una intersección física retenida dentro de esa asociación. Los sufijos ordinales se asignan en orden creciente de distancia sobre el perfil. Las etiquetas primer/intermedio/último se conservan en `etiqueta_pt`; pueden apuntar al mismo `id_pt`.

## Variables comunes

| Campo o familia | Significado | Unidad / convención |
| --- | --- | --- |
| `distancia_m`, `s_m` | Distancia acumulada del perfil exportado | m, creciente aguas arriba |
| `chi_m`, `chi_b_m`, `chi_pt_m` | Coordenada χ con θ=0,40 y A₀=10⁶ m² | m, origen declarado por perfil |
| `elevacion_m`, `z_m`, `z_kp_m` | Cota guardada del DEM | m; datum vertical absoluto pendiente |
| `area_m2` | Área aportante en el nodo | m² |
| `x_utm`, `y_utm`, sufijos `_m` | Coordenadas del cauce o selección derivada | m; WGS84 / UTM 18S, EPSG:32718 |
| `latitud_sur` | Magnitud positiva de latitud sur | grados |
| `H_dem_m`, `H_publicada_m` | Altura usada para buscar Pₜ y contraste publicado | m; H no es z(B) |
| `edad_mis_ka`, `error_edad_ka` | Edad asignada al nivel y rango publicado | ka |
| `K_G_m02_anio` | Erodabilidad efectiva por propagación | m⁰·²/año |
| `K_P_m02_anio` | Erodabilidad efectiva por pendiente y tasa costera | m⁰·²/año |
| `K_*_inferior_m02_anio`, `K_*_superior_m02_anio` | Coeficientes usados en el tiempo | m⁰·²/año; superior vacío con K uniforme |
| `tau_*_referencia_ka` | Tiempo de la versión final, conservado para cotejo | ka |
| `tau_*_ka` | Tiempo recalculado | ka |
| `chi_origen_m`, `z_origen_m`, `s_origen_m` | Inicio del tramo temporal c | m, mismas coordenadas que el perfil |
| `corte_ajuste_s_min_m` | Límite inferior del ajuste de revisión | m; puede diferir de c |
| `ksn_m08`, `ksn_B_PT_m08` | Empinamiento normalizado del ajuste indicado | m⁰·⁸ |
| `E_mm_anio` | 1.000 K_G ksn | mm/año, incisión modelada |
| `U_KP_mm_anio` | (z_KP − z_c) / τ_c | mm/año, alzamiento medio condicionado |
| `cota_absoluta_sobre_tau_mm_anio` | z_KP / τ_c | mm/año, contraste de referencia vertical |
| `R2`, `R2_B_PT` | Coeficiente de determinación del ajuste | Adimensional; no prueba equilibrio |
| `RMSE_geometrico_ka` | Dispersión de residuos B–Pₜ respecto de T | ka; no intervalo de confianza |

Los escenarios se identifican por la combinación río, MIS, método y modelo, no solo por el nombre del KP. Las ventanas conservan longitud solicitada, longitud efectiva y cobertura. Una ventana parcial no se presenta como completa.
