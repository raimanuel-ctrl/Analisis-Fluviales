"""Reproduce the final thesis calculations from the distributed derived inputs.

No DEM routing, spatial terrace matching or automatic KP acceptance is performed.
"""
from pathlib import Path
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / 'datos'
THETA, A0_M2 = 0.40, 1_000_000.0
FACTOR = A0_M2 ** THETA
SOUTH = [10, 12, 13, 15, 16, 18]


def leer(nombre):
    return pd.read_csv(DATA / nombre, float_precision='round_trip')


def perfil(rid):
    return pd.read_csv(DATA / 'perfiles' / f'R{int(rid):02d}_perfil.csv', float_precision='round_trip')


def calibraciones_centrales():
    """Count distinct physical PTs once; average K before computing time."""
    cal = leer('calibraciones_k_por_pt.csv')
    physical = cal.drop_duplicates(['rid', 'mis', 'id_pt'])
    return physical.groupby(['rid', 'mis'], as_index=False).agg(
        K_G_m02_anio=('K_G_m02_anio', 'mean'),
        K_P_m02_anio=('K_P_m02_anio', 'mean'),
        n_pt_fisicos=('id_pt', 'nunique'))


def tiempo(chi, origen, k_inferior, k_superior=None, frontera=None):
    """ka; chi and origin in the same exported profile coordinate system."""
    values = np.asarray(chi, dtype=float)
    if not np.isfinite(k_inferior) or k_inferior <= 0:
        raise ValueError('K inferior debe ser positivo y finito.')
    if k_superior is None or not np.isfinite(k_superior):
        return (values - origen) / (k_inferior * FACTOR * 1000)
    if k_superior <= 0 or frontera is None or not np.isfinite(frontera):
        raise ValueError('El modelo compuesto requiere K superior positivo y frontera finita.')
    def primitive(x):
        return (np.minimum(x, frontera) / k_inferior + np.maximum(x - frontera, 0) / k_superior) / (FACTOR * 1000)
    return primitive(values) - primitive(origen)


def tiempos_finales():
    """Recompute saved choices, with C2 only for Reloca."""
    t = leer('escenarios_tiempo.csv').copy()
    border = float(leer('reloca_fronteras.csv').query("frontera == 'C2'").chi_m.iloc[0])
    for method in ['G', 'P']:
        t[f'tau_{method}_ka'] = [float(tiempo(r.chi_kp_m, r.chi_origen_m,
            getattr(r, f'K_{method}_inferior_m02_anio'), getattr(r, f'K_{method}_superior_m02_anio'),
            border if r.rid == 6 else None)) for r in t.itertuples()]
    return t


def alzamiento_tramo(t=None):
    """Conditional mean-uplift upper bounds; metres/ka = millimetres/year."""
    t = tiempos_finales() if t is None else t
    output = []
    for r in t.itertuples():
        for method in ['G', 'P']:
            tau = getattr(r, f'tau_{method}_ka')
            output.append(dict(rid=r.rid, id_kp=r.id_kp, mis=r.mis, metodo=method,
                estado_kp=r.estado_kp, latitud_sur=r.latitud_sur,
                z_kp_m=r.z_kp_m, z_origen_m=r.z_origen_m,
                desnivel_m=r.z_kp_m-r.z_origen_m, tau_ka=tau,
                U_KP_mm_anio=(r.z_kp_m-r.z_origen_m)/tau,
                cota_absoluta_sobre_tau_mm_anio=r.z_kp_m/tau))
    return pd.DataFrame(output)


def resumen_sur(t=None):
    t = tiempos_finales() if t is None else t
    output = []
    for mis in ['5c', '5e']:
        for method in ['G', 'P']:
            a = t[t.rid.isin(SOUTH) & t.mis.eq(mis)][f'tau_{method}_ka']
            output.append(dict(mis=mis, metodo=method, n=len(a), minimo_ka=a.min(),
                maximo_ka=a.max(), promedio_ka=a.mean(), mediana_ka=a.median()))
    return pd.DataFrame(output)


def calibrar_reloca():
    """Recompute the individual positive C1-C4 lower-domain K values and residuals."""
    cal = leer('calibraciones_k_por_pt.csv').query('rid == 6').drop_duplicates('id_pt')
    means = calibraciones_centrales().query('rid == 9').set_index('mis')
    kp = leer('catalogo_kp.csv').query('rid == 6').iloc[0]
    origin = leer('origenes_temporales.csv').query('rid == 6').chi_origen_m.iloc[0]
    output = []
    for c in leer('reloca_fronteras.csv').itertuples():
        for mis in ['5c', '5e']:
            upper = means.loc[mis, 'K_G_m02_anio']
            ks = []
            for p in cal.itertuples():
                lower_delta = max(min(p.chi_pt_m, c.chi_m)-p.chi_b_m, 0.0)
                upper_delta = max(p.chi_pt_m-max(p.chi_b_m, c.chi_m), 0.0)
                remaining = p.edad_mis_ka*1000 - upper_delta/(upper*FACTOR)
                ks.append(lower_delta/(FACTOR*remaining) if remaining > 0 and lower_delta > 0 else np.nan)
            valid = bool(np.isfinite(ks).all() and np.all(np.asarray(ks) > 0))
            lower = float(np.mean(ks)) if valid else np.nan
            residuals = [float(tiempo(p.chi_pt_m, p.chi_b_m, lower, upper, c.chi_m))-p.edad_mis_ka for p in cal.itertuples()] if valid else []
            output.append(dict(frontera=c.frontera, mis_donante=mis, admisible=valid,
                K_inferior_m02_anio=lower, K_superior_m02_anio=upper,
                tau_kp_ka=float(tiempo(kp.chi_m, origin, lower, upper, c.chi_m)) if valid else np.nan,
                RMSE_geometrico_ka=float(np.sqrt(np.mean(np.square(residuals)))) if valid else np.nan))
    return pd.DataFrame(output)


def incision_ventanas():
    """The fits are saved inputs. Only the conversion K * ksn is repeated here."""
    d = leer('ajustes_incision_ventanas.csv').copy()
    d['E_mm_anio'] = 1000*d.K_G_m02_anio*d.ksn_m08
    return d


def verificar():
    t = tiempos_finales()
    for method in ['G', 'P']:
        np.testing.assert_allclose(t[f'tau_{method}_ka'], t[f'tau_{method}_referencia_ka'], rtol=1e-11, atol=1e-8)
    kp = leer('catalogo_kp.csv')
    assert len(kp) == 17 and kp.id_kp.is_unique
    assert set(kp.query("estado_kp == 'potencial'").rid) == {11, 14, 17}
    assert 7 not in set(kp.rid)
    for r in kp.itertuples():
        p = perfil(r.rid)
        np.testing.assert_allclose(np.interp(r.s_m, p.distancia_m, p.elevacion_m), r.z_m, atol=1e-5, rtol=0)
    means = calibraciones_centrales()
    assert int(means.n_pt_fisicos.sum()) == 18
    # Temporal calibration closes for individual PTs, not generally for their mean K.
    cal = leer('calibraciones_k_por_pt.csv')
    np.testing.assert_allclose((cal.chi_pt_m-cal.chi_b_m)/(cal.K_G_m02_anio*FACTOR*1000),cal.edad_mis_ka,atol=1e-8)
    for r in means.itertuples():
        use = t.query('rid == @r.rid and mis == @r.mis')
        if r.rid in [6, 7] or use.empty:
            continue
        np.testing.assert_allclose(use.K_G_inferior_m02_anio, r.K_G_m02_anio, rtol=1e-11)
        np.testing.assert_allclose(use.K_P_inferior_m02_anio, r.K_P_m02_anio, rtol=1e-11)
    c2 = calibrar_reloca().query("frontera == 'C2'").set_index('mis_donante')
    for r in t.query('rid == 6').itertuples():
        np.testing.assert_allclose(c2.loc[r.mis,'tau_kp_ka'], r.tau_G_ka, rtol=1e-11)
    d = incision_ventanas()
    np.testing.assert_allclose(d.E_mm_anio, d.E_referencia_mm_anio, atol=1e-10)
    s = resumen_sur(t).query("mis == '5e' and metodo == 'G'").iloc[0]
    np.testing.assert_allclose([s.minimo_ka,s.maximo_ka,s.promedio_ka], [459.98591851363926,1185.0913621200295,793.8155006626856],atol=1e-7)
    return dict(catalogo_kp=17,kp_potenciales=3,rio_sin_kp='R07',casos_calibracion=len(means),
        pt_fisicos=18,escenarios_tiempo=len(t),ajustes_incision=len(d),
        tiempos_coinciden_con_referencia=True,C2_recalculado=True,
        alcance='Reproduccion numerica desde derivados. No se reejecuta extraccion D8 desde el DEM.')


def figuras(out, t=None):
    import matplotlib.pyplot as plt
    t = tiempos_finales() if t is None else t
    kp = leer('catalogo_kp.csv')
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42}):
        fig, axs = plt.subplots(2,1,figsize=(10,7),layout='constrained')
        for rid in range(1,19):
            p = perfil(rid);ax=axs[0 if rid <= 9 else 1]
            ax.plot(p.chi_m/1000,p.elevacion_m,lw=.9,label=f'R{rid:02d}')
            for k in kp[kp.rid.eq(rid)].itertuples():
                ax.scatter(k.chi_m/1000,k.z_m,marker='*',s=75,facecolors='white' if k.estado_kp=='potencial' else '#262626',edgecolors='#262626',zorder=5)
        for ax in axs:
            ax.set(xlabel='χ desde el origen del perfil exportado (km)',ylabel='Elevación DEM (m)');ax.grid(alpha=.15);ax.legend(ncol=5,fontsize=8)
        _guardar(fig,out,'01_perfiles_chi_y_kp')
        fig, ax = plt.subplots(figsize=(10,4),layout='constrained')
        for mis,mark in [('5c','s'),('5e','o')]:
            q=t[t.mis.eq(mis)]
            for method,colour,offset in [('G','#126080',-.06),('P','#b55b0c',.06)]:
                for pot in [False,True]:
                    a=q[q.estado_kp.eq('potencial')==pot]
                    ax.scatter(a.rid+offset,a[f'tau_{method}_ka'],marker=mark,facecolors='white' if pot else colour,edgecolors=colour,label=f'K_{method}, MIS {mis}'+(' · potencial' if pot else ''))
        ax.set(xlabel='Río (orden norte a sur)',ylabel='τ del tramo adoptado (ka)',xticks=range(1,19));ax.legend(ncol=2,fontsize=8);ax.grid(alpha=.15)
        _guardar(fig,out,'02_tiempos_propagacion')
        u=alzamiento_tramo(t).query("metodo == 'G' and mis == '5e'")
        fig, ax = plt.subplots(figsize=(10,4),layout='constrained')
        for pot in [False,True]:
            a=u[u.estado_kp.eq('potencial')==pot]
            ax.scatter(a.latitud_sur,a.U_KP_mm_anio,facecolors='white' if pot else '#126080',edgecolors='#126080',label='KP potencial' if pot else 'KP reconocido')
            for r in a.itertuples():ax.annotate(f'R{r.rid:02d}',(r.latitud_sur,r.U_KP_mm_anio),xytext=(4,5),textcoords='offset points',fontsize=8)
        ax.set(xlabel='Latitud (°S) · norte a sur',ylabel='U_KP = Δz / τ (mm/año)',ylim=(0,None));ax.legend();ax.grid(alpha=.15)
        _guardar(fig,out,'03_alzamiento_medio_condicionado')


def _guardar(fig,out,name):
    import matplotlib.pyplot as plt
    for ext in ['png','svg']:
        fig.savefig(Path(out)/f'{name}.{ext}',dpi=250,facecolor='white')
    plt.close(fig)


def ejecutar(output=None, con_figuras=True):
    report=verificar()
    out=Path(output) if output else BASE/'resultados'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True,exist_ok=False)
    t=tiempos_finales()
    for name,df in [('calibraciones_k_centrales',calibraciones_centrales()),('tiempos_kp',t),
                    ('alzamiento_medio_tramo',alzamiento_tramo(t)),('resumen_seis_receptores',resumen_sur(t)),
                    ('reloca_fronteras_recalculadas',calibrar_reloca()),('incision_ventanas',incision_ventanas())]:
        df.to_csv(out/f'{name}.csv',index=False)
    if con_figuras:figuras(out,t)
    (out/'verificacion.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    return out,report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--sin-figuras',action='store_true')
    args=parser.parse_args()
    out,report=ejecutar(args.output,not args.sin_figuras)
    print(json.dumps(report,ensure_ascii=False,indent=2));print(out)
