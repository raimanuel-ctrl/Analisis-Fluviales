# External scientific data

The repository includes `notebooks/outlets_20261005.csv` for the 18 river outlets. Other scientific inputs/products are not distributed in this release.

| Branch | Required inputs |
| --- | --- |
| Topography | Original georeferenced DEM, outlets and regional JSON configuration. |
| Profiles | DEM and completed topography manifest/products. |
| Profile analysis | CSV profiles/manifest, manual KP CSV and optional candidate tables. |
| Terraces | Profiles, terrace workbook, sample-review workbook and selection workbook. |
| Faults/geology | Profiles, integrated polygons/faults, CHAF traces and optional terrace exports. |
| Detailed lithology | Profiles/manifest, manual KP CSV, polygons/faults, optional candidates/context tables. |
| Latest K/time | Profiles, manual KP CSV, `Ages_Final.xlsx`, terrace catalogue, declared outlet shapefile, dated lithology CSV/provenance JSON. |
| Regional interpretation | Valid upstream products and an explicitly populated scenario CSV. |

The regional JSON and detailed path blocks specify inputs. Shapefiles require associated files. Validate CRS and units; the current case uses EPSG:32718.

Profile geometry fields: `distancia_m`, `elevacion_m`, `area_m2`, `x_utm`, `y_utm`. Manifests record run IDs, exact paths and SHA-256 hashes; relocating files requires deliberate handling or fresh exports.

For distributed data, record authorship, source, citation, license, CRS, units and transformations. Explain how to obtain restricted inputs. Large data under `raw/` and `processed/` are ignored; a small redistributable example remains pending. Shareable manual selections and parameter decisions should be versioned, because automatic detection cannot reconstruct them.
