# Material computacional de la memoria

Versión actualizada el 10 de noviembre de 2026 a partir de la memoria final del 7 de noviembre. Esta carpeta reúne los notebooks, las selecciones manuales y los datos derivados que sostienen sus resultados. El material de investigación anterior se conserva en [archivo](../archivo/README.md).

## Empezar

Para reproducir los resultados finales, ejecutar [08 · Síntesis y tablas](notebooks/08_sintesis_y_tablas_memoria.ipynb). Para revisar el razonamiento, seguir el orden 01–08. El notebook 09 agrega la reproducción de los gráficos del informe. Cada notebook tiene una función y utiliza las mismas tablas de entrada.

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
| 09 | [Gráficos del informe](notebooks/09_graficos_del_informe.ipynb) | Reproducir las figuras por número, excluyendo los mapas y declarando las entradas externas. |

## Fórmulas y recorrido de los cálculos

Las ecuaciones siguientes describen la convención implementada en esta versión. Se distingue la geometría observada del perfil, las decisiones manuales de calibración y las cantidades calculadas a partir de ellas.

### 1. Referencias y unidades

| Símbolo | Significado | Unidad |
| --- | --- | --- |
| $s$ | Distancia acumulada sobre el perfil, creciente aguas arriba | m |
| $z$ | Elevación guardada del DEM | m |
| $A$ | Área aportante en cada nodo | m² |
| $A_0$ | Área de normalización, fijada en $10^6$ | m² |
| $\theta=m/n$ | Concavidad de referencia, fijada en 0,40 | Adimensional |
| $M_\chi$ | Pendiente del ajuste de $z$ frente a χ | Adimensional |
| $K$ | Erodabilidad efectiva con $m=0.40$ y $n=1$ | m⁰·²/año |
| $k_{sn}$ | Empinamiento normalizado con $\theta=0.40$ | m⁰·⁸ |
| $T$ | Edad MIS asignada a la referencia de terraza | ka, convertida a años para calibrar |
| $\tau_c$ | Tiempo modelado entre c y KP | ka |

**B** es el cruce reconstruido de la referencia geomorfológica con el río. **Pₜ** es una intersección del perfil actual con la elevación H de referencia. Puede haber varias. **c** fija el comienzo del recorrido temporal, mientras **KP** identifica el quiebre adoptado. La altura H de la terraza y la cota actual del cauce en B se registran por separado. La exageración vertical del visor afecta a la visualización, no a esas alturas.

### 2. Transformación χ y empinamiento

Con $s_0$ en el inicio del perfil exportado, la coordenada χ se define como:

$$
\chi(s)=\int_{s_0}^{s}\left(\frac{A_0}{A(u)}\right)^\theta\,du.
$$

La integración original utiliza trapecios sobre los nodos guardados:

$$
\chi_{i+1}-\chi_i=
\frac{s_{i+1}-s_i}{2}
\left[\left(\frac{A_0}{A_i}\right)^\theta+
\left(\frac{A_0}{A_{i+1}}\right)^\theta\right].
$$

La relación entre la pendiente del perfil transformado y el empinamiento normalizado es:

$$
z=\alpha+M_\chi\chi,
\qquad
k_{sn}=M_\chi A_0^\theta.
$$

Los ajustes lineales emplean intercepto libre. La pendiente depende del tramo o ventana seleccionados. Esta publicación conserva los ajustes de incisión guardados. El generador de sus gráficos reconstruye las rectas sobre los mismos intervalos. La posición de un KP puede haberse revisado mediante la intersección de dos rectas, manteniendo su elevación sobre el perfil observado.

### 3. Dos estimadores de K

**Propagación, K_G.** Para cada Pₜ físico se asigna la edad MIS al recorrido B–Pₜ y se despeja:

$$
K_{G,i}=\frac{\Delta\chi_{B-P_{t,i}}}{A_0^\theta T_{\mathrm{años}}}
=\frac{\chi(P_{t,i})-\chi(B)}{1000\,T_{\mathrm{ka}}\,A_0^\theta}.
$$

Los niveles utilizan $T=106\pm5$ ka para MIS 5c y $125\pm5$ ka para MIS 5e. La edad es una entrada de la calibración. Por eso, evaluar otra vez el mismo intervalo con su K individual recupera esa edad por construcción. Su significado físico requiere sostener la correlación de la terraza y la interpretación del recorrido.

**Pendiente y tasa costera, K_P.** La tasa de terraza se calcula como:

$$
U_{\mathrm{terraza},\mathrm{mm/año}}
=\frac{H_{\mathrm{publicada}}-z_{\mathrm{mar}}}{T_{\mathrm{ka}}}.
$$

Con alturas en m y edad en ka, m/ka es numéricamente equivalente a mm/año. La reconstrucción utilizada adopta $z_{\mathrm{mar}}=+5\pm5$ m para MIS 5e y $-20\pm5$ m para MIS 5c. Bajo el supuesto de incisión aproximadamente igual al alzamiento, se obtiene:

$$
K_P=\frac{U_{\mathrm{terraza},\mathrm{mm/año}}}{1000\,k_{sn,B-P_t}}.
$$

El factor 1.000 convierte la tasa a m/año. Como K_P incorpora la tasa costera, una comparación posterior con esa misma tasa comparte información de entrada.

**Valor central.** Para N intersecciones físicas distintas dentro de un mismo caso río/MIS:

$$
\overline K=\frac{1}{N}\sum_{i=1}^{N}K_i.
$$

Las etiquetas primer, intermedio y último pueden representar un mismo punto, por lo que `id_pt` evita contarlo varias veces. Se promedia K antes de calcular τ, manteniendo separados K_G, K_P y cada MIS. En general, el tiempo calculado con la media de K difiere de la media de los tiempos calculados por separado.

### 4. Tiempo de propagación con K uniforme

Para $n=1$, la velocidad de propagación del modelo en la coordenada χ es $K A_0^\theta$. El tiempo desde c hasta el KP es:

$$
\tau_{c,\mathrm{ka}}=
\frac{\chi(KP)-\chi(c)}{1000\,K\,A_0^\theta}.
$$

Los valores de `origenes_temporales.csv` definen c. Cambiar ese comienzo modifica el recorrido y su tiempo, mientras la calibración B–Pₜ se conserva y se controla por separado. Los segmentos anteriores a c pueden aparecer con coordenadas temporales negativas únicamente como contexto gráfico.

### 5. Reloca: integración con dos dominios de K

En Reloca, C2 separa los dominios inferior y superior. Su posición guardada es $\chi_C\approx2644.88$ m, a una elevación de 44,81 m. El dominio superior recibe $K_{sup}$ de Curanilahue, separando MIS 5c y 5e. La terraza local de Reloca mantiene 125 ka en ambos casos.

Para un Pₜ situado por encima de C2, el tiempo de calibración se descompone en:

$$
T_{\mathrm{años}}=
\frac{\chi_C-\chi_B}{K_{inf}A_0^\theta}
+\frac{\chi_{P_t}-\chi_C}{K_{sup}A_0^\theta}.
$$

De aquí se obtiene K inferior con el tiempo restante después de recorrer el dominio superior:

$$
K_{inf}=
\frac{\chi_C-\chi_B}
{A_0^\theta\left[T_{\mathrm{años}}-
\dfrac{\chi_{P_t}-\chi_C}{K_{sup}A_0^\theta}\right]}.
$$

La implementación recorta cada longitud al dominio efectivamente recorrido. Si Pₜ queda antes de C2, la longitud superior es cero y el recorrido inferior termina en Pₜ. La alternativa requiere longitudes y tiempo restante que entreguen un K inferior positivo. Esta regla también se aplica al evaluar C1–C4.

Para acumular el tiempo entre dos posiciones cualesquiera se define una primitiva y se resta su valor en c:

$$
F(\chi)=\frac{1}{1000\,A_0^\theta}
\left[\frac{\min(\chi,\chi_C)}{K_{inf}}+
\frac{\max(\chi-\chi_C,0)}{K_{sup}}\right],
\qquad
\tau_c=F(\chi_{KP})-F(\chi_c).
$$

La síntesis utiliza el promedio de los K inferiores admisibles por Pₜ físico. Los residuos de B–Pₜ respecto de 125 ka se expresan mediante:

$$
r_i=\tau(B,P_{t,i})-125,
\qquad
\mathrm{RMSE}_{geom}=\sqrt{\frac{1}{N}\sum_{i=1}^{N}r_i^2}.
$$

Ese RMSE describe dispersión geométrica en ka. C2 se mantiene como reconstrucción geológica preferida, mientras C3 y C4 presentan menor dispersión en este control.

### 6. Incisión y alzamiento medio del tramo

La ley de incisión utilizada y el balance de elevación son:

$$
E=K A^m S^n,
\qquad
\frac{\partial z}{\partial t}=U-E.
$$

Con $n=1$, $m=\theta=0.40$ y tasas en m/año, $E=Kk_{sn}$. Por tanto, la tasa de incisión de un ajuste guardado se presenta como:

$$
E_{\mathrm{mm/año}}=1000\,K_G\,k_{sn}.
$$

La equivalencia con alzamiento exige un balance suficientemente ajustado y un nivel de base representativo. Para el relieve del tramo se usa una estimación distinta, con el mismo inicio y final que el tiempo:

$$
U_{KP,\mathrm{mm/año}}=
\frac{z_{KP}-z_c}{\tau_{c,\mathrm{ka}}}.
$$

Su interpretación como cota superior condicionada supone una pendiente inicial mucho menor que la actual. Una contribución positiva de relieve heredado reduce la fracción atribuible al desplazamiento tectónico. El cociente alternativo $z_{KP}/\tau_c$ se conserva para comparar la referencia vertical.

Si K es uniforme, puede escribirse la identidad geométrica:

$$
U_{KP}=1000\,K_G\,A_0^\theta\frac{\Delta z}{\Delta\chi}
=1000\,K_G\,k_{sn,sec}.
$$

Aquí $k_{sn,sec}$ utiliza la secante entre los extremos del tramo. El ksn ajustado en una ventana próxima al KP puede ser diferente. En Reloca se mantiene la integración compuesta de dos K.

### 7. Ejemplo con un resultado guardado

Curanipe (R10) recibe el promedio de Curanilahue–MIS 5e:

$$
\overline K_G=1.1049607893\times10^{-5}\ \mathrm{m^{0.2}/año},
\qquad A_0^{0.40}\approx251.188643.
$$

Con $\chi_c=641.851734$ m y $\chi_{KP}=2860.742623$ m, la diferencia es 2218,890889 m. Sustituyendo los valores de precisión completa se obtiene:

$$
\tau_c\approx799.45\ \mathrm{ka}.
$$

Las cotas $z_c=22.105705$ m y $z_{KP}=373.469902$ m entregan un desnivel de 351,364197 m, de modo que:

$$
U_{KP}\approx\frac{351.364197}{799.445901}
=0.43951\ \mathrm{mm/año}.
$$

El ejemplo corresponde a los datos de `escenarios_tiempo.csv` y se recalcula en `src/analisis.py`. La elección principal del sur utiliza MIS 5e por su representación en el registro de terrazas considerado, mientras que MIS 5c permanece como contraste.

Las relaciones de χ y propagación se apoyan en [Perron y Royden (2013)](https://doi.org/10.1002/esp.3302) y [Royden y Perron (2013)](https://doi.org/10.1002/jgrf.20031). El modelo de incisión y su balance se desarrollan en [Whipple y Tucker (1999)](https://doi.org/10.1029/1999JB900120). Las edades y reconstrucciones costeras corresponden al uso de [Jara-Muñoz et al. (2015)](https://doi.org/10.1016/j.quascirev.2015.01.005) descrito en la memoria; la asignación de T al recorrido B–Pₜ conserva su supuesto específico de calibración.

## Herramientas interactivas

- [Explorador de KP](visores/Explorador_KP.html): archivo HTML autónomo con perfiles, mapa, tributarios de área ≥1 km² y regresiones exploratorias. La distribución abre con un borrador vacío. Los puntos adoptados se identifican mediante `KP-R##` y los clics nuevos permanecen exploratorios.
- [Visor de B en 3D](visores/Visor_B_3D.html): relieve y bandas N–S/O–E para los ocho B adoptados. Conserva las funciones de giro, dimensiones, exageración vertical y selección manual del visor utilizado en la investigación. La versión pública omite las ubicaciones originales del catálogo de terrazas.

Descargar los HTML para utilizarlos en el navegador. La vista de archivos de GitHub no ejecuta aplicaciones. Ambos funcionan sin servidor. Para regenerarlos desde esta carpeta:

```sh
python -m pip install -r requirements_visores.txt
python -m src.visores kp
python -m src.visores b3d
```

El visor 3D requiere el DEM de `data/dem/`, descargado mediante Git LFS. Su ruta se configura en [configuracion.json](configuracion.json). Las referencias de terrazas pueden añadirse desde un CSV autorizado localmente. El HTML generado con ese catálogo incorpora sus puntos y debe conservarse localmente.

El explorador calcula un K por medición. El modelo compuesto C2 se reproduce en los notebooks y no se aplica automáticamente en esa aplicación.

## Datos y tablas

[datos](datos/README.md) contiene los perfiles, el catálogo de KP, las selecciones B/Pₜ, las calibraciones y los escenarios con precisión numérica. [tablas_memoria](tablas_memoria/README.md) reúne trece tablas transcritas de la versión final, con su numeración y redondeo de presentación. El [diccionario de datos](docs/DICCIONARIO_DATOS.md) define los campos y las unidades. Las tablas redondeadas sirven para cotejar el texto. Los cálculos leen los datos de precisión completa.

La reproducción parte de datos derivados. La extracción D8 original desde el DEM y las asociaciones espaciales con terrazas requieren sus insumos y su entorno original, conservados o identificados en el archivo metodológico. La comprobación numérica de esta carpeta no constituye una nueva validación geomorfológica de las selecciones.

## Convenciones científicas

θ = 0,40; n = 1; A₀ = 10⁶ m²; K en m⁰·²/año. K_G designa propagación y K_P, pendiente y tasa costera. Se promedia K sobre Pₜ físicos distintos dentro de cada caso río/MIS, antes de calcular τ. El caso principal del sur transfiere K_G de Curanilahue–5e, con su cuenca completa. MIS 5c permanece como contraste.

Los tiempos parten de c, el comienzo del tramo temporal adoptado. B delimita la calibración con Pₜ. Estas referencias se guardan por separado. U_KP = (z_KP − z_c)/τ_c se interpreta como una cota superior condicionada del alzamiento medio bajo una pendiente inicial mucho menor que la actual. K_G ksn representa incisión modelada. Las edades y las tasas conservan los supuestos del modelo y el datum vertical absoluto del DEM sigue sin confirmación.

R01 conserva un KP sin K adoptado. R07 tiene calibraciones de K y carece de KP aceptado. R11, R14 y R17 son potenciales. La síntesis principal del sur incluye R10, R12, R13, R15, R16 y R18.

Ver [reproducibilidad](docs/METODOLOGIA_Y_REPRODUCIBILIDAD.md), [selección de contenidos](docs/SELECCION_DE_CONTENIDOS.md) y [cómo referenciar el repositorio en la memoria](docs/REPOSITORIO_EN_LA_MEMORIA.md).
