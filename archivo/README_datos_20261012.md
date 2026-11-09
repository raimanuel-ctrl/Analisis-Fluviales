# Published and external scientific data

Spatial datasets published on 2026-10-12:

- `dem/TDM_cut_UTM18.tif`: working DEM, stored with Git LFS. 512,582,580 bytes; 12 m grid; EPSG:32718; one float32 band; 14,077 × 12,815 cells; nodata=-99999.
- `outlets/outlets_20261005.csv`: 18 point-outlet records, identical to `notebooks/outlets_20261005.csv`.
- `outlets/rios_outlets.shp` and sidecars: the exact layer declared in `Calibracion_K_t`. It contains 18 main-channel LineStrings, not 18 points. Source export identifier: `20260929_144303_906995`.
- `geology/Geologia_unidades.shp` and sidecars: the supplied integrated map, 656 valid polygons/multipolygons in EPSG:32718. The homonymous QML and alternative Jorquera-color QML are included.
- Complete ZIP bundles are provided for both shapefiles. Keep `.shp`, `.shx`, `.dbf`, `.prj` and `.cpg` together. QGIS can load the homonymous style; the alternative dated QML can be loaded manually.

Files were copied without altering their bytes and checked against originals using SHA-256. `data_manifest_20261012.json` records file sizes/hashes; `metadata_20261012.json` records inspected spatial metadata. Geometry validity does not prove cartographic accuracy or geological correlation.

## Download and configure

```sh
git lfs install
git clone https://github.com/raimanuel-ctrl/Analisis-Fluviales.git
cd Analisis-Fluviales
git lfs pull
```

Git LFS is needed for the DEM. The repository's default ZIP may contain only its pointer; download the actual raster from its GitHub file page if using browser downloads. Shapefile ZIP bundles use normal Git storage.

For the regional workflow, set `dem` in `notebooks/configuracion_post_auditoria_20261005.json` to the downloaded DEM's absolute path and `geologia_integrada` to the downloaded `Geologia_unidades.shp`. Keep `base_dir` pointing to the intended input/product workspace. Detailed notebooks also require updating their `BASE`, `GEO_DIR`, `RUTA_GEO` and `RUTA_OUTLETS` path blocks where applicable. `rios_outlets` is a declared main-channel source, not a replacement for the point catalogue.

## Source limitations

`geology/LEEME.txt` and `equivalencias.csv` retain the compilation's source/correlation notes. The mosaicking rule prioritizes Jorquera > Chanco > Gajardo without establishing geological equivalence; map seams are not automatically contacts. Chanco classifications and some regional correlations remain preliminary. The source notes also describe fault/GPKG layers that are not included in this distribution. The code's MIT license does not confer rights over upstream datasets; no blanket data license is assigned here.

## Inputs still required

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

The table describes workflow requirements; the DEM, catalogue and geological-unit polygons above are now supplied. Exported hydrology/profiles, age/selection workbooks, KP CSVs and additional fault layers remain external. Regional JSON and detailed path blocks specify the remaining inputs.

Profile geometry fields: `distancia_m`, `elevacion_m`, `area_m2`, `x_utm`, `y_utm`. Manifests record run IDs, exact paths and SHA-256 hashes; relocating files requires deliberate handling or fresh exports.

For distributed data, record authorship, source, citation, license, CRS, units and transformations. Explain how to obtain restricted inputs. Large data under `raw/` and `processed/` are ignored; a small redistributable example remains pending. Shareable manual selections and parameter decisions should be versioned, because automatic detection cannot reconstruct them.
