# Technical context for scientific code review

Snapshot: 2026-10-12. This reference complements the root README with a reading sequence, terminology and criteria for scientific code review. It describes the implementation rather than newly executed results.

## Suggested reading order

1. `README.md`: objective, branches, parameters, equations and limits.
2. `ai_context/Calibracion_K_t.py`: latest sources/scenarios, crossing fits, K, outlet-to-KP times and K8 sensitivity.
3. `ai_context/Litologia_perfiles_comparaciones.py`: profile loading, manual rematching, geology and slope/ksn comparisons.
4. `notebooks/post_auditoria_20261005.py` plus its configuration JSON: regional provenance, QC and definitions.
5. `ai_context/Interpretaciones_K_deltaT_20261005.py`: gated regional interpretation, using a different time origin.
6. `docs/REVIEW_20261012.md`: review evidence and unresolved dependency.

Code-only files retain notebook cell boundaries and may include IPython magics; they are not standalone scripts. Notebook files are authoritative. The DEM (Git LFS), point-outlet catalogue, declared main-channel shapefile and integrated geological-unit map are now distributed under `data/`. Saved numerical outputs, workbooks and remaining scientific inputs are external. Consult the data README for paths, source limitations and hashes.

## Terminology and distinctions

- `perfil`: river profile, oriented from outlet upstream.
- `cruce`: river-profile elevation crossing, possibly one of several; not the spatial projection of a terrace.
- `terraza`: terrace marker; a correlated marker is not necessarily independently dated.
- `muestra`: dated sample/sample-rate source.
- `manual_reubicado`: detailed XY/Z-matched manual point, not automatically accepted by the regional hash policy.
- `completo`: export state, not final scientific validation.
- `pendiente`/`rechazado`: pending/rejected, not zero or physical absence.
- θ=0.45 is regional; θ=0.40 is recomputed in detailed notebooks using inherited geometry.
- ksn = Mχ A0^θ; K = (U/1000)/ksn for n=1. Units depend on θ.
- Regional Δt starts at the terrace crossing; detailed time starts at the loaded outlet. Check actual trimming externally.
- K8's K1/K2/K3 mean outlet–crossing, crossing–KP and outlet–KP fits. They use the same selected rate and equilibrium assumption.
- K8 refits each raw/smoothed curve independently but does not update earlier base K/time tables.

## Questions for scientific/code review

1. Are sample–terrace association, uplift-rate transfer, datum/sea-level treatment and equilibrium assumptions defensible?
2. Which crossings and bedrock intervals remain defensible after reviewing alluvium, geometry and DEM limits?
3. Do mapped lithological changes prevent constant K transfer toward the KP?
4. What uncertainty remains outside crossing/smoothing sensitivity?
5. Are normalization, time origin, regression conventions and units consistent with the exact implementation used for any external comparison?
6. Which downstream products require regeneration when inputs or decisions change?

The algebra of the K core can be inspected here. Full parity with another researcher's workflow requires its exact module, filters, regression, geometry and uncertainty treatment; similarity of equations alone is insufficient.

Keep code inspection, locally recorded exports, reproducibility and scientific interpretation distinct. Do not infer numerical results from stripped notebooks. No hydrology, overlay, detector or K/time calculation was executed for this publication.
