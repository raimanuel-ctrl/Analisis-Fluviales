# Fluvial profiles and coastal landscape evolution in central Chile

Python workflows for river morphometry, coastal terraces and conditional stream-power erodibility calibration. Geology thesis research by Raimundo Sánchez Mosqueira, Universidad de Chile.

## Current version

The **22 October 2026** update is in [notebooks/actual_20261022](notebooks/actual_20261022/README.md). It contains seven ordered notebooks, supporting modules and an explicit derived-data snapshot for the current thesis synthesis.

Start with [07 Synthesis of current results](notebooks/actual_20261022/Notebooks_principales_20261019/07_Sintesis_resultados_20261022.ipynb) to reproduce mean K, propagation times, Reloca C2 calibration residuals and modeled incision rates immediately below knickpoints. It includes its derived inputs and produces static PNG/SVG figures and CSV tables; no interactive application or unpublished external comparison is needed.

The current synthesis has been numerically checked against saved thesis tables. The six upstream research notebooks have been checked for schema and syntax but were **not rerun** for this release. Full reproduction from the DEM still needs additional external inputs and a separately validated environment. These levels of verification are not interchangeable.

## Scientific scope

Eighteen coastal channels between 35°S and 36.5°S are compared using longitudinal profiles, χ, normalized steepness and mapped geology/faults. Six channels provide eight river/MIS calibration cases. The current analysis uses θ = 0.40, n = 1 and A₀ = 10⁶ m²; dimensional K is in m⁰·²/year.

K_G denotes propagation calibration and K_P denotes slope/rate calibration. Means count distinct physical terrace-elevation crossings once and keep MIS 5c and 5e separate. Reloca uses the conditional C2 two-domain model, with residual checks at the original calibration points. Drainage-area and time-origin alternatives remain explicit. Modeled propagation times are not direct tectonic dates, and modeled incision is equated to uplift only under local equilibrium.

The current results do not confirm a unique 0.8–1 Ma uplift pulse. Calibration method and donor MIS change the absolute time scale. Likewise, the lack of an unequivocal recent fluvial response near Pelluhue and Carranza does not establish fault inactivity.

## Repository layout

| Location | Contents |
| --- | --- |
| [Current workflow](notebooks/actual_20261022/README.md) | Seven notebooks, supporting modules, dependency lists, provenance and current derived inputs. |
| [Current numerical synthesis](notebooks/actual_20261022/sintesis_20261022.py) | Standalone equivalent of notebook 07. Creates dated local outputs without overwriting inputs. |
| [Spatial data](data/README.md) | Working DEM through Git LFS, 18-outlet catalogue, main-channel layer and integrated geological compilation from the earlier release. |
| [Verification notes](docs/REVIEW_20261022.md) | Static checks, executed synthesis, exclusions and remaining reproduction limits. |
| [Earlier workflow](docs/SNAPSHOT_20261012.md) | Documentation of the 12 October release; its ten notebooks and code-only reading snapshots remain available for provenance. Not the current source of K/τ results. |

The earlier regional workflow used θ = 0.45 in some stages. Do not mix its K units, χ coordinates or outputs with the current θ = 0.40 calculations. Detailed differences and the required execution order are documented in each snapshot.

## Reuse and publication boundaries

Install the synthesis dependencies, open Jupyter in the current directory and run notebook 07 in order, or run its standalone module. See the [current README](notebooks/actual_20261022/README.md) for commands and the complete workflow. Large rasters require `git lfs install` and `git lfs pull`; a source ZIP is not assumed to contain LFS raster bytes.

Source notebooks are published without saved outputs, execution counts or personal filesystem paths. Private documents, correspondence, collaborators' original notebooks and unpublished external-comparison values are not redistributed. Omitted inputs are identified rather than replaced with invented examples. Original code is covered by the [MIT license](LICENSE); scientific datasets and third-party software retain their own conditions and attribution.
