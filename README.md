# Fluvial profiles, coastal terraces and conditional erodibility — central Chile

Python workflows for fluvial morphometry, coastal-terrace analysis and conditional calibration of stream-power erodibility in central Chile. Developed for thesis research in Geology at Universidad de Chile.

**Author:** Raimundo Sanchez Mosqueira. **Current version:** research snapshot dated 12 October 2026. **Status:** active research development.

The current implementation includes the revised regional workflow and subsequent lithological and calibration analyses. The `20261005` filenames identify the workflow generation rather than the date of every subsequent edit. Earlier versions are preserved in Git history.

## Scientific objective and current direction

The project investigates relationships among river-profile morphology, knickpoints, coastal terrace markers, mapped geology and landscape evolution in central Chile. The regional workflow builds profiles for 18 outlets. Detailed calibration focuses on six rivers and evaluates the effects of reference elevation, sample uplift rate, crossing choice and smoothing on conditional K estimates and modeled propagation times.

The central question is whether a downstream calibration can be transferred toward a manually selected knickpoint, given lithological changes, alluvial reaches, structural context and terrace-correlation uncertainty. Spatial coincidence and a good χ–elevation fit alone do not establish that transfer or a tectonic origin.

**This is ongoing research, not a validated tectonic chronology.** The repository includes code, the working DEM, the outlet catalogue, the declared river shapefile and the integrated geological-unit map. Exported profiles, additional fault layers and age workbooks remain external. This publication review did not execute the notebooks or recalculate results.

## Published spatial data

| Dataset | Location | Contents |
| --- | --- | --- |
| Working DEM | [TDM_cut_UTM18.tif](data/dem/TDM_cut_UTM18.tif) | Original bytes, 512,582,580 bytes, 12 m pixels, EPSG:32718; stored through Git LFS. |
| Point-outlet catalogue | [outlets_20261005.csv](data/outlets/outlets_20261005.csv) | 18 outlet coordinates used by regional topography; identical to the copy beside the notebooks. |
| Declared `rios_outlets` shapefile | [Download complete ZIP](data/outlets/rios_outlets_20261012.zip) | The exact layer referenced by detailed calibration: 18 main-channel LineStrings, despite its name; not a point-outlet layer. |
| Integrated geological units | [Download complete ZIP](data/geology/Geologia_unidades_20261012.zip) | 656 polygons/multipolygons from the Chanco–Jorquera–Gajardo compilation, all shapefile sidecars, QGIS styles and source notes. |

With Git LFS installed, `git clone` and `git lfs pull` retrieve the actual DEM. Without it, the repository may contain a small pointer rather than raster bytes. See [GitHub's Git LFS documentation](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-git-large-file-storage). For browser downloads, use the DEM page's raw/download control. A default repository ZIP should not be assumed to include LFS data.

See [data instructions](data/README.md), [metadata](data/metadata_20261012.json) and [file hashes](data/data_manifest_20261012.json). Update notebook/configuration paths to the downloaded datasets; published paths are not automatically adopted by the research notebooks. The integrated map retains preliminary classifications, scale differences and map-seam limitations documented in `data/geology/LEEME.txt`; broad lithological groups do not establish identical erodibility. Scientific datasets do not automatically inherit the code's MIT license; source conditions and attribution still apply.

## Technical documentation

The detailed analysis is implemented in [Calibracion_K_t.ipynb](notebooks/Calibracion_K_t.ipynb), which contains the 30-scenario K/time analysis and K8 crossing/smoothing sensitivity. [Litologia_perfiles_comparaciones.ipynb](notebooks/Litologia_perfiles_comparaciones.ipynb) provides lithological assignments and regional profile comparisons.

[AI_CONTEXT.md](AI_CONTEXT.md) documents terminology, methodological distinctions and review criteria. The [ai_context/](ai_context/) directory provides code-only notebook snapshots with cell boundaries for text-based inspection and automated review. These files are reading aids rather than standalone scripts; notebooks remain authoritative. The shared Python module and JSON configuration document the regional infrastructure.

## Workflow and files

| Branch | Notebook | Role |
| --- | --- | --- |
| Regional 1 | [AnalisisTopografico_20261005](notebooks/AnalisisTopografico_20261005.ipynb) | DEM/CRS checks, D8 hydrology, area, 18 basins, QC and manifests. |
| Regional 2 | [Perfiles_chi_ksn_20261005](notebooks/Perfiles_chi_ksn_20261005.ipynb) | Main-channel extraction, longitudinal profiles, trapezoidal χ, local ksn and window sensitivity. |
| Regional 3 | [AnalisisPerfiles_20261005](notebooks/AnalisisPerfiles_20261005.ipynb) | Regional fits, local candidates, optional TopoToolbox/external methods, manual selections and rejection/coincidence tables. |
| Regional 4 | [AnalisisTerrazas_20261005](notebooks/AnalisisTerrazas_20261005.ipynb) | Terrace associations, projections, elevation crossings and vertical-uncertainty alternatives; Excel/CSV selection reconciliation. |
| Spatial branch | [Perfiles_y_fallas_20261005](notebooks/Perfiles_y_fallas_20261005.ipynb) | Integrated geology, contacts, river–fault intersections/proximities; optional CHAF transects. Source updated on 12 October. |
| Regional interpretation | [Interpretaciones_K_deltaT_20261005](notebooks/Interpretaciones_K_deltaT_20261005.ipynb) | Explicit decision table and gated conditional K/Δt scenarios. Empty decisions remain pending. |
| Detailed lithology | [Litologia_perfiles_comparaciones](notebooks/Litologia_perfiles_comparaciones.ipynb) | χ recalculated at θ=0.40; manual-point relocation, mapped units, slope/ksn and concavity sensitivity. |
| Detailed calibration | [Calibracion_K_t](notebooks/Calibracion_K_t.ipynb) | Six rivers, 30 sample-rate/reference-elevation combinations, K, outlet-to-KP times, diagnostic fits, lithology and K8 sensitivity. |
| Supplementary mapping | [Litologia_mapa_knickpoints_20261005](notebooks/Litologia_mapa_knickpoints_20261005.ipynb) | Earlier focused mapping experiment, including candidates and proposed alluvial intervals. |
| Alternative calibration | [Analisis_K_rios_seleccionados_20261005](notebooks/Analisis_K_rios_seleccionados_20261005.ipynb) | Earlier scenarios comparing apparent relief/age and published-rate calibration, with explicit alluvial exclusions; separate from the latest 30-scenario calculation. |

Keep the regional notebooks beside [post_auditoria_20261005.py](notebooks/post_auditoria_20261005.py), [configuration](notebooks/configuracion_post_auditoria_20261005.json) and [outlets](notebooks/outlets_20261005.csv). Run regional stages 1–4 in order. The spatial branch can start after profiles; optional terrace overlays require stage 4. Interpretation needs the relevant exports and explicit decisions. Detailed notebooks load exported products in their own fresh kernel.

## Geometry, parameters and provenance

The published DEM header was read on 12 October: 12 m pixels, EPSG:32718, 14,077 rows × 12,815 columns, one float32 band and nodata=-99999. No hydrology or numerical profile analysis was rerun. The updated code checks actual resolution rather than using the old fixed `30m` label.

Distance increases upstream from the loaded outlet. Elevation/distance use metres and area uses m². Main-channel extraction follows the upstream candidate with greatest contributing area; original DEM elevations define profiles, while conditioned products support routing.

| Parameter | Regional workflow | Detailed lithology/calibration |
| --- | --- | --- |
| θ | 0.45 | 0.40; χ recomputed from geometry |
| A0 | 1,000,000 m² | 1,000,000 m² |
| n in interpretation | 1 | 1 |
| Channel threshold | 1,000,000 m² | Inherited profiles; no fresh channel extraction |
| Initial local slope window | 500 m; 250/500/1000 m sensitivity | Separate slope and smoothing sensitivities |

Local ksn uses a distance–elevation slope window; calibration uses a χ–elevation segment fit. Non-positive slopes are not converted with an absolute value. These quantities, θ values and window settings must be distinguished before comparing results.

The regional run is `20261005_post_auditoria_01`. Products are dated under `Productos/post_auditoria/<run_id>/<stage>`, indexed by `manifest_<run_id>.json`. Manifests contain paths, SHA-256 hashes, configuration and export states. `completo` means the stage recorded an export, not that its interpretation is final or its dependencies remain current.

Regional manual-selection loading validates hashes and indices. Detailed notebooks rematch manual KP coordinates/elevations with XY and Z tolerances of 0.10 m, keeping original/current hashes and rejected points. These acceptance policies differ.

## Equations and time origins

```text
χ(s) = integral from outlet to s of (A0/A(s'))^θ ds'  [trapezoidal integration]
z = b + Mχ χ                                        [free-intercept regression]
ksn_segment = Mχ A0^θ
K = (U_m_ka/1000) / ksn_segment                      [n=1, U assumed equal to erosion]
t_ka = Δχ / (K A0^θ) / 1000                         [constant K, stationary area]
```

For θ=0.40, dimensional K has units m^0.2/year; for θ=0.45, m^0.1/year. Mχ and dimensional ksn are different quantities. Do not apply another A0^θ correction to K after using this formula.

**Time origins differ.** Regional interpretation calibrates between an explicit bedrock start and the terrace crossing, then computes Δχ from that crossing to the KP. `Calibracion_K_t` computes time from the loaded outlet (`s=0`) to the manual KP. Its labels say trimmed outlet, but the notebook does not itself trim/reconstruct profiles. Input origin and alluvial coverage must be independently checked. These times are not directly interchangeable or automatically tectonic ages.

## Six-river scenarios and K8

| ID | River | Sample-rate alternatives in latest calibration |
| --- | --- | --- |
| 3 | Estero Quebrada Honda | CONT 11 |
| 4 | Estero Malhueco | CONT 11 |
| 7 | Estero de Chanco | PUHE 10, PUHE 11; exploratory case |
| 8 | Río Rahue | PUHE 10, PUHE 11 |
| 9 | Río Curanilahue | PUHE 10, PUHE 11 |
| 10 | Río Curanipe | PUHE 10, PUHE 11; correlations need review |

Each sample-rate alternative is combined with two configured terrace elevations and the sample's own elevation: 30 scenario definitions. `Ages_Final.xlsx` supplies ages, sample elevations and `Urate (m/ka)`; a separate catalogue supplies terrace elevations. These combinations are sensitivity scenarios, not proof that each terrace was dated with that sample. OD values are separately transcribed in code, not read from `Ages_Final`. Primary age/rate provenance was not reverified in this publication review.

The base calculation defaults to the first elevation crossing unless `CRUCES_ELEGIDOS` overrides it. It fits outlet–crossing with at least five original nodes and an interpolated endpoint, records alternatives and rejects unsuitable fits. Five nodes are a numerical threshold, not a quality guarantee.

In K8, **K1 = outlet–crossing, K2 = crossing–KP, K3 = outlet–KP**. These are three fitting intervals using the selected rate, not independent physical measurements. All available crossings are examined, with the first marked principal in each model. Multiple manual KPs require explicit selection.

K8 supports raw elevation, Savitzky–Golay or moving-average smoothing and optional joint display. Each curve gets its own crossings, fits and K1/K2/K3. Smoothing preserves distance, area, χ and KP position; elevations are filtered on a regular distance grid, with effective odd-node window and border treatment recorded. **K8 does not replace the earlier base K table or recompute its travel-time table with smoothed elevations.**

## Geological context and unresolved work

Mapped units/faults provide context. Bilateral river-bank queries can resolve narrow mapping gaps using a configurable 20 m search limit; this is an operational setting, not cartographic accuracy. Detailed controls flag rivers 3 and 4 as lacking confirmed map coverage. Seams, ambiguous overlaps and proximity remain distinct from demonstrated contacts or tectonic controls.

This release passed JSON/Jupyter-schema and Python syntax checks for ten notebooks and the shared module. Public scientific cell sources match the originals after declared personal-path replacement. See [review notes](docs/REVIEW_20261012.md) and [publication hashes](docs/publication_provenance_20261012.json).

Read-only local manifest checks matched checked files in topography, profiles, profile analysis, terraces and faults; four large entries were not rehashed. **Interpretation has a changed upstream terrace-manifest dependency and must be regenerated.** This is provenance evidence, not a numerical rerun.

Remaining work:

- Regenerate interpretation after reviewing upstream terrace exports and scenario decisions.
- Justify sample–terrace correlations, rate transfer, datum/sea-level assumptions and alluvial/bedrock limits.
- Check loaded profile origins and geological compatibility along the extrapolation reach.
- Update the explicit lithology date if new products are adopted; K7 currently uses `20261006_132435_105931`.
- Validate optional method adapters and calculations in the intended environment.
- Document exact package sources/versions and a redistributable executable example.

## Installation and reuse

Start Jupyter in `notebooks/` and review paths before running cells. Personal paths are replaced by examples under `C:\datos`. The regional JSON controls common paths; detailed notebooks contain extra path blocks. External manifests store machine-specific paths/hashes; moving them requires deliberate provenance handling or fresh exports.

Dependencies include NumPy, pandas, SciPy, Matplotlib, Rasterio, GeoPandas, Shapely, IPython/Jupyter, an Excel engine such as openpyxl and TopoAnalysis. TopoToolbox is optional for candidate comparison. An exact validated environment is pending; no untested lockfile is supplied.

See [data requirements](data/README.md). Original code contributions are distributed under the [MIT license](LICENSE). External software and datasets retain their respective licenses and attribution. TopoAnalysis and TopoToolbox are third-party dependencies. Bugs, methodological questions and proposed improvements can be submitted through repository issues.
