# Perfiles fluviales y evolución del relieve costero de Chile central

Material computacional de la memoria de título en Geología de Raimundo Sánchez Mosqueira, Universidad de Chile. Incluye los notebooks, instrumentos de selección y datos derivados utilizados para comparar los perfiles de 18 cauces entre 35°S y 36,5°S.

**Versión de entrega: 2026.11.10.** Los resultados corresponden a la memoria final del 7 de noviembre. La carpeta [memoria](memoria/README.md) es la única entrada al flujo vigente.

## Accesos principales

| Material | Acceso |
| --- | --- |
| Notebooks ordenados | [Ruta de trabajo 01–09](memoria/README.md#empezar) |
| Fórmulas y supuestos | [De χ y ksn a K, τ y alzamiento](memoria/README.md#fórmulas-y-recorrido-de-los-cálculos) |
| Gráficos del informe, sin mapas | [Cobertura por figura y comando de reproducción](memoria/docs/GRAFICOS_DEL_INFORME.md) |
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

Para abrir una aplicación, descargar el HTML y abrirlo en un navegador. GitHub muestra el archivo como código. La distribución funciona localmente sin enviar datos a un servidor.

## Reproducir la síntesis

```sh
git clone https://github.com/raimanuel-ctrl/Analisis-Fluviales.git
cd Analisis-Fluviales/memoria
python -m pip install -r requirements.txt
python -m src.analisis
```

Cada ejecución crea una carpeta nueva con seis tablas numéricas, tres figuras de control en PNG/SVG y la verificación de resultados. La síntesis utiliza datos derivados incluidos. El DEM es necesario para regenerar el relieve 3D y para las etapas topográficas originales. Permanece en `data/dem/` mediante Git LFS.

Para generar los gráficos del informe, conservando su numeración y excluyendo los mapas:

```sh
python -m src.graficos_informe
```

El inventario identifica las 56 figuras del documento. Ocho son mapas que quedan fuera de este generador. De la figura mixta 14 se reproduce únicamente el perfil. Las entradas de las comparaciones externas se aportan localmente, como detalla la [revisión de cobertura](memoria/docs/GRAFICOS_DEL_INFORME.md).

## Cómo se conectan los cálculos

El análisis parte de los perfiles del cauce, donde se conocen elevación, distancia y área aportante. La transformación χ permite expresar la geometría sobre una coordenada que incorpora el cambio del área aguas arriba:

$$
\chi(s)=\int_{s_0}^{s}\left(\frac{A_0}{A(u)}\right)^{\theta}\,du,
\qquad k_{sn}=M_\chi A_0^\theta.
$$

Aquí, $s_0$ es el inicio del perfil exportado y $M_\chi$ la pendiente ajustada de elevación frente a χ. Se utilizan $\theta=0.40$, $A_0=10^6\ \mathrm{m^2}$ y $n=1$. B y Pₜ delimitan la calibración con la terraza. El punto c es el inicio del recorrido temporal y KP, el knickpoint adoptado. Son referencias con funciones diferentes.

Asignando al recorrido B–Pₜ la edad de referencia $T$, el estimador de propagación y el tiempo para K uniforme se expresan como:

$$
K_G=\frac{\chi(P_t)-\chi(B)}{1000\,T_{\mathrm{ka}}\,A_0^\theta},
\qquad
\tau_{c,\mathrm{ka}}=\frac{\chi(KP)-\chi(c)}{1000\,K\,A_0^\theta}.
$$

El factor 1.000 convierte años y ka. Los promedios de K cuentan una vez cada Pₜ físico y mantienen separados los MIS. Reloca utiliza la suma de tiempos de dos dominios, descrita con detalle en el [README del flujo](memoria/README.md#5-reloca-integración-con-dos-dominios-de-k).

La incisión modelada y la estimación del alzamiento medio del tramo conservan sus propias definiciones:

$$
E_{\mathrm{mm/año}}=1000\,K_G\,k_{sn},
\qquad
U_{KP,\mathrm{mm/año}}=\frac{z_{KP}-z_c}{\tau_{c,\mathrm{ka}}}.
$$

La primera usa la pendiente del ajuste elegido. La segunda utiliza el desnivel y el tiempo entre los mismos extremos, bajo el supuesto geomorfológico indicado a continuación. El [desarrollo completo](memoria/README.md#fórmulas-y-recorrido-de-los-cálculos) incluye K_P, tasas de terraza, el modelo C2 y un ejemplo numérico de Curanipe.

## Alcance científico

Se utilizan θ = 0,40, n = 1 y A₀ = 10⁶ m². Las ocho asociaciones río/MIS reúnen 18 Pₜ físicos distintos. El promedio central de K cuenta cada punto una vez y conserva separados MIS 5c y 5e. Reloca utiliza C2 con dos dominios efectivos. Los receptores del sur reciben el K de la cuenca completa de Curanilahue, con MIS 5e como escenario principal.

Los tiempos de propagación se calculan desde los comienzos de tramo adoptados. U_KP = (z_KP − z_c)/τ_c se interpreta como una cota superior condicionada del alzamiento medio bajo una pendiente inicial mucho menor que la actual. Las tasas de incisión, las tasas costeras y U_KP conservan sus definiciones y duraciones de integración. La reproducción del cálculo no demuestra por sí sola una edad tectónica, un pulso regional único ni la inactividad de una falla.

R01 conserva KP sin K adoptado. R07 tiene K y carece de KP aceptado. R11, R14 y R17 permanecen potenciales. Los seis receptores principales del sur son R10, R12, R13, R15, R16 y R18.

## Organización y conservación

- [memoria/](memoria/README.md): versión vigente, datos derivados, herramientas y tablas finales.
- [data/](data/README.md): DEM, drenaje y geología publicados previamente, con sus condiciones de origen.
- [archivo/](archivo/README.md): versiones de octubre y antecedentes metodológicos, conservados para trazabilidad.

Las versiones anteriores quedan identificadas como históricas para evitar mezclar K, geometrías y orígenes temporales. Los originales de trabajo se conservan fuera de esta reorganización. Las [notas de selección](memoria/docs/SELECCION_DE_CONTENIDOS.md) explican qué material se mantuvo activo.

El código propio está cubierto por [MIT](LICENSE). Las fuentes científicas conservan su autoría y condiciones. Los catálogos espaciales privados de terrazas, Ages_Final, correspondencia, manuscritos y valores no publicados de colaboradores quedan fuera de la distribución. La reproducción integral desde el DEM requiere insumos adicionales y el entorno de las etapas originales.

## Citar

Sánchez Mosqueira, R. (2026). *Analisis-Fluviales: material computacional de la memoria de geomorfología fluvial* (versión 2026.11.10) [Código y datos derivados]. GitHub. https://github.com/raimanuel-ctrl/Analisis-Fluviales

Los metadatos están en [CITATION.cff](CITATION.cff). Para la entrega de la memoria se debe conservar también el identificador del commit utilizado.
