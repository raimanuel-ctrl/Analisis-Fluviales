"""Swath del polígono guardado y esquema 3D a escala real.

Importar desde Jupyter. Este archivo no ejecuta cálculos al importarse.
No modifica el DEM, shapefiles ni notebooks.
"""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from functools import lru_cache
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import rasterio
from affine import Affine
from matplotlib.colors import LightSource
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds
from scipy.ndimage import binary_dilation
from shapely.geometry import LineString, Point


def _partes_lineales(geom):
    if geom.is_empty:
        return
    if geom.geom_type == "LineString":
        yield geom
    elif hasattr(geom, "geoms"):
        for parte in geom.geoms:
            yield from _partes_lineales(parte)


def validar_ksn(perfiles_ksn, ventana_m=None, theta=None):
    """Comprobar los parámetros efectivamente usados, no solo el encabezado."""
    if not perfiles_ksn:
        raise ValueError("Primero ejecutar la celda de cálculo de ksn.")
    primero = next(iter(perfiles_ksn.values()))
    if ventana_m is None:
        ventana_m = primero.attrs.get("ventana_m", np.nan)
    if theta is None:
        theta = primero.attrs.get("theta_ksn_mapa", np.nan)
    ventana_m, theta = float(ventana_m), float(theta)
    if not np.isfinite(ventana_m) or ventana_m <= 0:
        raise ValueError("La ventana debe ser positiva y finita.")
    if not np.isfinite(theta) or theta < 0:
        raise ValueError("Theta debe ser finito y no negativo.")
    for rid, p in perfiles_ksn.items():
        actual = p.attrs.get("ventana_m", np.nan)
        if not np.isfinite(actual) or not np.isclose(actual, ventana_m):
            raise ValueError(
                f"Río {rid}: datos calculados con {actual:g} m, "
                f"pero el encabezado solicita {ventana_m:g} m. "
                "Ejecuta nuevamente la celda de cálculo de ksn."
            )
        if not np.isclose(p.attrs.get("theta_ksn_mapa", np.nan), theta):
            raise ValueError(f"Río {rid}: recalcular ksn con theta={theta:g}.")
    return ventana_m, theta


def preparar_swath(ruta_poligono, ruta_dem, paso_m=250.0,
                   lado_max_3d=900, bloque=512):
    """Estadísticas a resolución nativa; malla reducida solo para dibujar."""
    if not np.isfinite(paso_m) or paso_m <= 0 or lado_max_3d < 10:
        raise ValueError("Revisar paso_m y lado_max_3d.")
    capa = gpd.read_file(ruta_poligono)
    if capa.crs is None or len(capa) != 1:
        raise ValueError("Se requiere un único polígono guardado con CRS.")
    pol = capa.to_crs(32718).geometry.iloc[0]
    if pol.geom_type != "Polygon" or not pol.is_valid or pol.is_empty:
        raise ValueError("El polígono no es válido.")

    rect = np.asarray(pol.minimum_rotated_rectangle.exterior.coords)[:-1, :2]
    lados = np.roll(rect, -1, axis=0) - rect
    u = lados[np.argmax(np.linalg.norm(lados, axis=1))]
    u = u / np.linalg.norm(u)
    if np.dot(u, [1, 1]) < 0:
        u = -u
    v = np.array([-u[1], u[0]])
    vertices = np.asarray(pol.exterior.coords)[:, :2]
    origen = vertices[np.argmin(vertices @ u)]
    s_max = float(((vertices - origen) @ u).max())
    t_min, t_max = ((vertices - origen) @ v).min(), ((vertices - origen) @ v).max()
    nb = int(np.ceil(s_max / paso_m))
    cuenta = np.zeros(nb, dtype=np.int64)
    suma = np.zeros(nb, dtype=np.float64)
    minimo = np.full(nb, np.inf)
    maximo = np.full(nb, -np.inf)

    with rasterio.open(ruta_dem) as dem:
        if dem.crs is None or dem.crs.to_epsg() != 32718:
            raise ValueError("El DEM debe estar en EPSG:32718.")
        if dem.transform.b != 0 or dem.transform.d != 0:
            raise ValueError("Se requiere una grilla sin rotación.")
        xmin, ymin, xmax, ymax = pol.bounds
        completa = (xmin >= dem.bounds.left and xmax <= dem.bounds.right
                    and ymin >= dem.bounds.bottom and ymax <= dem.bounds.top)
        w = from_bounds(xmin, ymin, xmax, ymax, transform=dem.transform)
        c0, r0 = max(0, int(np.floor(w.col_off))), max(0, int(np.floor(w.row_off)))
        c1 = min(dem.width, int(np.ceil(w.col_off + w.width)))
        r1 = min(dem.height, int(np.ceil(w.row_off + w.height)))
        if c1 <= c0 or r1 <= r0:
            raise ValueError("El polígono no intersecta el DEM.")
        win = Window(c0, r0, c1 - c0, r1 - r0)
        total = int(np.ceil(win.width / bloque) * np.ceil(win.height / bloque))
        procesados = 0
        print(f"Swath sobre DEM original: {dem.res[0]:g} m; {total} bloques.")
        for fila in range(r0, r1, bloque):
            for col in range(c0, c1, bloque):
                tw = Window(col, fila, min(bloque, c1-col), min(bloque, r1-fila))
                tf = dem.window_transform(tw)
                z = dem.read(1, window=tw, masked=True, out_dtype="float32")
                dentro = geometry_mask([pol.__geo_interface__],
                                      out_shape=z.shape, transform=tf, invert=True)
                valido = dentro & ~np.ma.getmaskarray(z) & np.isfinite(z.data)
                if valido.any():
                    xx = tf.c + (np.arange(z.shape[1]) + .5) * tf.a
                    yy = tf.f + (np.arange(z.shape[0]) + .5) * tf.e
                    ss = (xx[None, :] - origen[0])*u[0] + (yy[:, None] - origen[1])*u[1]
                    indices = np.clip(np.floor(ss[valido]/paso_m).astype(int), 0, nb-1)
                    valores = z.data[valido].astype(float)
                    cuenta += np.bincount(indices, minlength=nb)
                    suma += np.bincount(indices, weights=valores, minlength=nb)
                    np.minimum.at(minimo, indices, valores)
                    np.maximum.at(maximo, indices, valores)
                procesados += 1
                if procesados % 25 == 0 or procesados == total:
                    print(f"Bloques procesados: {procesados}/{total}")

        # Malla para la superficie 3D. No interviene en las estadísticas.
        factor = min(1.0, lado_max_3d / max(win.width, win.height))
        ny, nx = max(2, int(win.height*factor)), max(2, int(win.width*factor))
        zg = dem.read(1, window=win, out_shape=(ny, nx),
                      masked=True, out_dtype="float32").filled(np.nan)
        tg = dem.window_transform(win) * Affine.scale(win.width/nx, win.height/ny)
        mascara = geometry_mask([pol.__geo_interface__], out_shape=zg.shape,
                               transform=tg, invert=True)
        zg[~mascara] = np.nan
        xg = tg.c + (np.arange(nx)+.5)*tg.a
        yg = tg.f + (np.arange(ny)+.5)*tg.e
        resolucion = tuple(dem.res)

    valido = cuenta > 0
    media = np.divide(suma, cuenta, out=np.full(nb, np.nan), where=valido)
    minimo[~valido], maximo[~valido] = np.nan, np.nan
    desde = np.arange(nb)*paso_m
    hasta = np.minimum(desde+paso_m, s_max)
    tabla = pd.DataFrame({
        "desde_m": desde, "hasta_m": hasta, "centro_m": (desde+hasta)/2,
        "z_min_m": minimo, "z_media_m": media, "z_max_m": maximo,
        "n_pixeles": cuenta,
    })
    if not np.isfinite(zg).any():
        raise ValueError("La malla visual no contiene elevaciones válidas.")
    fill = np.where(np.isfinite(zg), zg, np.nanmedian(zg))
    shade = LightSource(azdeg=315, altdeg=45).hillshade(
        fill, vert_exag=1, dx=abs(tg.a), dy=abs(tg.e))
    shade[binary_dilation(~np.isfinite(zg), iterations=1)] = np.nan
    return {
        "poligono": pol, "origen": origen, "u": u, "v": v,
        "azimut": float(np.degrees(np.arctan2(u[0], u[1])) % 360),
        "s_max_m": s_max, "t_min_m": float(t_min), "t_max_m": float(t_max),
        "swath": tabla, "x_m": xg, "y_m": yg, "z_m": zg, "hillshade": shade,
        "ruta_dem": str(Path(ruta_dem)), "resolucion_dem_m": resolucion,
        "poligono_dentro_dem": bool(completa),
        "paso_swath_m": paso_m,
    }


def figura_swath_3d(modelo, perfiles_ksn, segmentos_ksn, manuales,
                    fallas_chaf, fallas_maule=None, mostrar_resumen=True,
                    nombres=None, ventana_m=None, theta=None):
    """Relieve y trazas reales. No dibuja planos de falla supuestos."""
    ventana_m, theta = validar_ksn(perfiles_ksn, ventana_m, theta)
    unidades_ksn = f"m^{2*theta:g}"
    nombres = {} if nombres is None else nombres
    pol = modelo["poligono"]
    xmin, ymin, xmax, ymax = pol.bounds
    x0, y0 = xmin, ymin
    fig = go.Figure()
    fig.add_trace(go.Surface(
        x=(modelo["x_m"]-x0)/1000, y=(modelo["y_m"]-y0)/1000,
        z=modelo["z_m"]/1000, surfacecolor=modelo["hillshade"],
        colorscale=[[0, "#555555"], [1, "#f3f3f3"]],
        cmin=0, cmax=1, showscale=False, opacity=.88,
        lighting=dict(ambient=1, diffuse=0, specular=0),
        name="DEM dentro de banda", connectgaps=False,
        hovertemplate="Elevación: %{z:.3f} km s.n.m.<extra>DEM</extra>",
    ))

    todos_k = np.concatenate([d["ksn_040"] for d in segmentos_ksn.values()])
    buenos = todos_k[np.isfinite(todos_k) & (todos_k > 0)]
    if not len(buenos):
        raise ValueError("No hay ksn válidos.")
    cmin, cmax = np.percentile(buenos, [2, 98])
    if cmax <= cmin:
        cmin, cmax = cmin/2, cmax*2
    lmin, lmax = np.log10([cmin, cmax])
    ticks = np.unique(np.r_[cmin, [k for k in [1,2,5,10,20,50,100,200,500,1000]
                                  if cmin < k < cmax], cmax])
    ticktext = [f"{k:.1f}" for k in ticks]
    ticktext[0], ticktext[-1] = "≤"+ticktext[0], "≥"+ticktext[-1]

    xyz_validos, colores, meta, xyz_invalidos, etiquetas, filas = [], [], [], [], [], []
    for rid in sorted(segmentos_ksn, key=int):
        p, d = perfiles_ksn[rid], segmentos_ksn[rid]
        if len(d["xy"]) != len(p)-1:
            raise ValueError(f"Río {rid}: segmentos y perfil incompatibles.")
        seg = d["xy"]
        elegibles = ((seg[:,:,0].max(axis=1) >= xmin) & (seg[:,:,0].min(axis=1) <= xmax)
                     & (seg[:,:,1].max(axis=1) >= ymin) & (seg[:,:,1].min(axis=1) <= ymax))
        z = p.elevacion_m.to_numpy(float)
        longitud, n = 0., 0
        puntos_rio = []
        for j in np.flatnonzero(elegibles):
            xy_original = seg[j]
            inter = LineString(xy_original).intersection(pol)
            delta = xy_original[1]-xy_original[0]
            den = float(np.dot(delta, delta))
            if den <= 0 or not np.isfinite(z[j:j+2]).all():
                continue
            for parte in _partes_lineales(inter):
                xy = np.asarray(parte.coords)[:, :2]
                fraccion = np.clip(((xy-xy_original[0]) @ delta)/den, 0, 1)
                zz = z[j] + fraccion*(z[j+1]-z[j])
                xyz = np.column_stack([(xy[:,0]-x0)/1000, (xy[:,1]-y0)/1000, zz/1000])
                puntos_rio.extend(xyz)
                longitud += parte.length
                n += 1
                kval = float(d["ksn_040"][j])
                if np.isfinite(kval) and kval > 0:
                    xyz_validos.extend(xyz.tolist()+[[np.nan]*3])
                    colores.extend([np.log10(kval)]*(len(xyz)+1))
                    meta.extend([[rid, nombres.get(rid, rid), kval]]*(len(xyz)+1))
                else:
                    xyz_invalidos.extend(xyz.tolist()+[[np.nan]*3])
        if puntos_rio:
            etiquetas.append((rid, puntos_rio[len(puntos_rio)//2]))
            filas.append({"ID": rid, "rio": nombres.get(rid, rid),
                          "longitud_en_banda_km": longitud/1000,
                          "segmentos_visibles": n})

    if xyz_invalidos:
        a = np.asarray(xyz_invalidos)
        fig.add_trace(go.Scatter3d(x=a[:,0], y=a[:,1], z=a[:,2], mode="lines",
                                  line=dict(color="#888888", width=4),
                                  name="Ríos: sin ksn válido", hoverinfo="skip"))
    if xyz_validos:
        a = np.asarray(xyz_validos)
        fig.add_trace(go.Scatter3d(
            x=a[:,0], y=a[:,1], z=a[:,2], mode="lines", connectgaps=False,
            customdata=np.asarray(meta, dtype=object),
            line=dict(color=colores, width=6, colorscale="Turbo",
                      cmin=lmin, cmax=lmax, showscale=True,
                      colorbar=dict(title=f"ksn ({unidades_ksn})<br>color logarítmico",
                                    tickvals=np.log10(ticks), ticktext=ticktext, len=.65)),
            name=f"Ríos: ksn local {ventana_m:g} m",
            hovertemplate=("%{customdata[1]} (ID %{customdata[0]})<br>"
                           f"ksn: %{{customdata[2]:.2f}} {unidades_ksn}<br>"
                           f"Ventana: {ventana_m:g} m<br>"
                           "Elevación: %{z:.3f} km<extra></extra>"),
        ))
    if etiquetas:
        a = np.asarray([q for rid, q in etiquetas])
        fig.add_trace(go.Scatter3d(x=a[:,0], y=a[:,1], z=a[:,2], mode="text",
                                  text=[rid for rid, q in etiquetas],
                                  textfont=dict(size=12, color="black"),
                                  name="ID de ríos", showlegend=False, hoverinfo="skip"))

    # Muestreo de trazas sobre el DEM nativo por bloques, sin rellenar huecos.
    with rasterio.open(modelo["ruta_dem"]) as dem:
        @lru_cache(maxsize=64)
        def leer_bloque(br, bc):
            r, c = br*256, bc*256
            return dem.read(1, window=Window(c, r, min(256, dem.width-c),
                                             min(256, dem.height-r)),
                            masked=True, out_dtype="float32")

        def elevaciones(xy):
            resultado = np.full(len(xy), np.nan)
            rr, cc = rasterio.transform.rowcol(dem.transform, xy[:,0], xy[:,1])
            rr, cc = np.asarray(rr), np.asarray(cc)
            validos = (rr >= 0) & (rr < dem.height) & (cc >= 0) & (cc < dem.width)
            inds = np.flatnonzero(validos)
            for j in inds:
                val = leer_bloque(int(rr[j]//256), int(cc[j]//256))[rr[j]%256, cc[j]%256]
                if not np.ma.is_masked(val) and np.isfinite(val):
                    resultado[j] = float(val)
            return resultado

        fuentes = [(fallas_chaf, False)]
        if fallas_maule is not None:
            fuentes.append((fallas_maule, True))
        for capa, contexto in fuentes:
            for _, falla in capa.to_crs(32718).iterrows():
                for parte in _partes_lineales(falla.geometry.intersection(pol)):
                    dist = np.linspace(0, parte.length,
                                       max(2, int(np.ceil(parte.length/150))+1))
                    xy = np.array([(parte.interpolate(s).x, parte.interpolate(s).y) for s in dist])
                    zz = elevaciones(xy)
                    nombre = "Maule: contexto" if contexto else str(falla.get("falla", "CHAF"))
                    if not contexto:
                        nombre += " | " + str(falla.get("id_traza", falla.get("f_id", "")))
                    fig.add_trace(go.Scatter3d(
                        x=(xy[:,0]-x0)/1000, y=(xy[:,1]-y0)/1000, z=zz/1000,
                        mode="lines", connectgaps=False,
                        line=dict(color="#555555" if contexto else "#81246a",
                                  width=5, dash="dash" if contexto else "solid"),
                        name=nombre, hovertemplate=nombre+"<extra></extra>",
                    ))
        dist = np.linspace(0, pol.exterior.length, max(2, int(pol.exterior.length/200)+1))
        xy = np.array([(pol.exterior.interpolate(s).x, pol.exterior.interpolate(s).y) for s in dist])
        zz = elevaciones(xy)
        fig.add_trace(go.Scatter3d(x=(xy[:,0]-x0)/1000, y=(xy[:,1]-y0)/1000,
                                  z=zz/1000, mode="lines", line=dict(color="black", width=3),
                                  name="Límite de banda", connectgaps=False))
        leer_bloque.cache_clear()

    kp = manuales.loc[[pol.covers(Point(r.x_utm, r.y_utm)) for r in manuales.itertuples()]]
    if not kp.empty:
        fig.add_trace(go.Scatter3d(
            x=(kp.x_utm-x0)/1000, y=(kp.y_utm-y0)/1000, z=kp.elevacion_m/1000,
            mode="markers", marker=dict(size=5, color="black"), name="KP manuales",
            text=kp.id_kp, hoverinfo="text",
        ))

    # Perfil estadístico al costado SE de la banda: no es una traza del terreno.
    if mostrar_resumen:
        tabla = modelo["swath"]
        margen = max(1000., .06*(modelo["t_max_m"]-modelo["t_min_m"]))
        xy = (modelo["origen"] + tabla.centro_m.to_numpy()[:,None]*modelo["u"]
              + (modelo["t_min_m"]-margen)*modelo["v"])
        for columna, etiqueta, color in [
            ("z_min_m", "Swath: mínimo", "#2a9d8f"),
            ("z_media_m", "Swath: media", "#111111"),
            ("z_max_m", "Swath: máximo", "#e76f51"),
        ]:
            fig.add_trace(go.Scatter3d(
                x=(xy[:,0]-x0)/1000, y=(xy[:,1]-y0)/1000,
                z=tabla[columna]/1000, mode="lines", line=dict(color=color, width=5),
                customdata=np.column_stack([tabla.centro_m/1000, tabla[columna], tabla.n_pixeles]),
                name=etiqueta+" (resumen lateral)", connectgaps=False,
                hovertemplate="Distancia SW→NE: %{customdata[0]:.2f} km<br>"
                              "Elevación: %{customdata[1]:.1f} m<br>"
                              "Píxeles: %{customdata[2]:.0f}<extra>Swath agregado</extra>",
            ))

    camara = dict(eye=dict(x=-1.8, y=-2.3, z=.5),
                  projection=dict(type="orthographic"), up=dict(x=0, y=0, z=1))
    fig.update_layout(
        title=f"Banda del graben | ksn local {ventana_m:g} m | theta={theta:g} | escala real 1:1",
        height=850, margin=dict(l=0, r=0, b=0, t=90),
        scene=dict(aspectmode="data", camera=camara,
                   xaxis_title="Este relativo (km)", yaxis_title="Norte relativo (km)",
                   zaxis_title="Elevación (km s.n.m.)"),
        legend=dict(font=dict(size=10)),
        updatemenus=[dict(type="buttons", direction="right", x=0, y=1.05, buttons=[
            dict(label="Sur–SO → NE", method="relayout", args=[{"scene.camera": camara}]),
            dict(label="Desde el sur", method="relayout", args=[{"scene.camera": dict(
                eye=dict(x=0,y=-2.8,z=.5), projection=dict(type="orthographic"),
                up=dict(x=0,y=0,z=1))}]),
            dict(label="Planta", method="relayout", args=[{"scene.camera": dict(
                eye=dict(x=0,y=0,z=3), projection=dict(type="orthographic"),
                up=dict(x=0,y=1,z=0))}]),
        ])],
    )
    return fig, pd.DataFrame(filas)


def mapa_ksn_interactivo(perfiles_ksn, segmentos_ksn, manuales, nombres,
                         ruta_dem, fallas_chaf, fallas_maule=None,
                         ventana_m=None, theta=None, rios=None,
                         max_pixeles_hillshade=1800, ancho_cauces=5,
                         n_colores=40):
    """Mapa con hillshade y zoom; parámetros tomados del cálculo de ksn."""
    import io
    import base64
    from PIL import Image
    from plotly.colors import sample_colorscale

    VENTANA_MAPA_KSN_M, THETA_MAPA_KSN = validar_ksn(
        perfiles_ksn, ventana_m=ventana_m, theta=theta
    )
    ksni_unidades = f"m^{2 * THETA_MAPA_KSN:g}"
    gf_ksn_perfiles, gf_ksn_segmentos = perfiles_ksn, segmentos_ksn
    gf_manuales, gf_nombres = manuales, nombres
    if fallas_chaf.crs is None:
        raise ValueError("Las trazas CHAF requieren un CRS.")
    gf_trazas_mapa = fallas_chaf.to_crs(32718).copy()
    if "id_traza" not in gf_trazas_mapa:
        gf_trazas_mapa["id_traza"] = [
            f"CHAF:{i}" for i in gf_trazas_mapa.index
        ]
    if fallas_maule is None:
        gf_maule_mapa = gf_trazas_mapa.iloc[:0].copy()
    else:
        if fallas_maule.crs is None:
            raise ValueError("Las trazas Maule requieren un CRS.")
        gf_maule_mapa = fallas_maule.to_crs(32718)
    RIOS_INTERACTIVOS_KSN = (
        list(sorted(segmentos_ksn, key=int)) if rios is None
        else [str(i) for i in rios]
    )
    if not RIOS_INTERACTIVOS_KSN:
        raise ValueError("Selecciona al menos un río.")
    for rid in RIOS_INTERACTIVOS_KSN:
        if rid not in perfiles_ksn or rid not in segmentos_ksn:
            raise ValueError(f"No hay datos calculados para el río {rid}.")
    MAX_PIXELES_HILLSHADE = int(max_pixeles_hillshade)
    ANCHO_CAUCES_KSN = float(ancho_cauces)
    N_COLORES_KSN = int(n_colores)
    if MAX_PIXELES_HILLSHADE < 2 or N_COLORES_KSN < 2 or ANCHO_CAUCES_KSN <= 0:
        raise ValueError("Revisar resolución, colores o ancho de los cauces.")
    ksni_xy_total = np.vstack([
        gf_ksn_perfiles[i][["x_utm", "y_utm"]].to_numpy(float)
        for i in RIOS_INTERACTIVOS_KSN
    ])
    ksni_min = ksni_xy_total.min(axis=0) - 3000
    ksni_max = ksni_xy_total.max(axis=0) + 3000

    with rasterio.open(Path(ruta_dem)) as ksni_dem:
        if ksni_dem.crs is None or ksni_dem.crs.to_epsg() != 32718:
            raise ValueError("El DEM debe estar en EPSG:32718.")

        if ksni_dem.transform.b != 0 or ksni_dem.transform.d != 0:
            raise ValueError("Esta vista requiere una grilla sin rotación.")

        ksni_bounds = (
            max(ksni_min[0], ksni_dem.bounds.left),
            max(ksni_min[1], ksni_dem.bounds.bottom),
            min(ksni_max[0], ksni_dem.bounds.right),
            min(ksni_max[1], ksni_dem.bounds.top)
        )

        if (
            ksni_bounds[0] >= ksni_bounds[2]
            or ksni_bounds[1] >= ksni_bounds[3]
        ):
            raise ValueError("Los ríos no intersectan el DEM.")

        ksni_win = from_bounds(
            *ksni_bounds, transform=ksni_dem.transform
        )
        ksni_win = ksni_win.intersection(
            Window(0, 0, ksni_dem.width, ksni_dem.height)
        )

        ksni_ratio = min(1.0, MAX_PIXELES_HILLSHADE / max(ksni_win.width, ksni_win.height))
        ksni_ny = max(2, int(ksni_win.height * ksni_ratio))
        ksni_nx = max(2, int(ksni_win.width * ksni_ratio))

        ksni_z = ksni_dem.read(
            1,
            window=ksni_win,
            out_shape=(ksni_ny, ksni_nx),
            masked=True,
            out_dtype="float32"
        ).filled(np.nan)

        (
            ksni_left, ksni_bottom,
            ksni_right, ksni_top
        ) = ksni_dem.window_bounds(ksni_win)

    ksni_valid = np.isfinite(ksni_z)
    if not ksni_valid.any():
        raise ValueError("No hay elevaciones válidas para el hillshade.")

    # El relleno se usa únicamente para calcular la iluminación.
    # Los huecos y sus bordes se vuelven transparentes en la imagen.
    ksni_fill = np.where(
        ksni_valid, ksni_z, np.nanmedian(ksni_z)
    )
    ksni_shade = LightSource(
        azdeg=315, altdeg=45
    ).hillshade(
        ksni_fill,
        vert_exag=1,
        dx=(ksni_right - ksni_left) / ksni_nx,
        dy=(ksni_top - ksni_bottom) / ksni_ny
    )

    ksni_alpha = ~binary_dilation(
        ~ksni_valid, iterations=1
    )
    ksni_rgba = np.zeros(
        (ksni_ny, ksni_nx, 4), dtype=np.uint8
    )
    ksni_rgba[:, :, :3] = (
        255 * np.clip(ksni_shade, 0, 1)
    ).astype(np.uint8)[:, :, None]
    ksni_rgba[:, :, 3] = ksni_alpha.astype(np.uint8) * 255

    ksni_buffer = io.BytesIO()
    Image.fromarray(ksni_rgba).save(
        ksni_buffer, format="PNG"
    )
    ksni_fondo = (
        "data:image/png;base64,"
        + base64.b64encode(
            ksni_buffer.getvalue()
        ).decode()
    )

    print("Hillshade preparado.")
    import plotly.graph_objects as go
    from plotly.colors import sample_colorscale

    # ---------- Segmentos y datos para consultar con el cursor ----------
    ksni_segments, ksni_values, ksni_meta = [], [], []

    for ksni_id in RIOS_INTERACTIVOS_KSN:
        ksni_d = gf_ksn_segmentos[ksni_id]

        ksni_segments.append(ksni_d["xy"])
        ksni_values.append(ksni_d["ksn_040"])
        ksni_meta.append(np.column_stack([
            np.full(
                len(ksni_d["ksn_040"]),
                gf_nombres[ksni_id],
                dtype=object
            ),
            ksni_d["ksn_040"],
            ksni_d["s0_m"],
            ksni_d["s1_m"],
            np.where(
                ksni_d["ventana_parcial"],
                "parcial", "completa"
            )
        ]))

    ksni_segments = np.concatenate(ksni_segments)
    ksni_values = np.concatenate(ksni_values)
    ksni_meta = np.concatenate(ksni_meta)

    ksni_ok = np.isfinite(ksni_values) & (ksni_values > 0)
    if not ksni_ok.any():
        raise ValueError("No hay ksn positivos disponibles.")

    # Escala logarítmica común para mejorar el contraste.
    # P2 y P98 son límites visuales; los valores calculados se conservan.
    ksni_cmin, ksni_cmax = np.percentile(
        ksni_values[ksni_ok], [2, 98]
    )
    if ksni_cmax <= ksni_cmin:
        ksni_cmin, ksni_cmax = ksni_cmin / 2, ksni_cmax * 2

    ksni_lmin, ksni_lmax = np.log10([
        ksni_cmin, ksni_cmax
    ])

    ksni_bin = np.full(len(ksni_values), -1, dtype=int)
    ksni_bin[ksni_ok] = np.clip(
        np.floor(
            (np.log10(ksni_values[ksni_ok]) - ksni_lmin)
            / (ksni_lmax - ksni_lmin)
            * N_COLORES_KSN
        ).astype(int),
        0, N_COLORES_KSN - 1
    )

    ksni_colors = sample_colorscale(
        "Turbo",
        (np.arange(N_COLORES_KSN) + 0.5) / N_COLORES_KSN
    )
    ksni_fig = go.Figure()

    # ---------- Cauce completo ----------
    ksni_gray = np.concatenate([
        np.vstack([
            gf_ksn_perfiles[i][["x_utm", "y_utm"]].to_numpy(float),
            [[np.nan, np.nan]]
        ])
        for i in RIOS_INTERACTIVOS_KSN
    ])

    ksni_fig.add_trace(go.Scattergl(
        x=ksni_gray[:, 0] / 1000,
        y=ksni_gray[:, 1] / 1000,
        mode="lines",
        line=dict(color="#888888", width=3),
        name="Cauce completo; gris = sin ksn válido",
        hoverinfo="skip"
    ))

    # ---------- Segmentos coloreados ----------
    for ksni_b, ksni_color in enumerate(ksni_colors):
        ksni_sel = ksni_bin == ksni_b
        if not ksni_sel.any():
            continue

        # Separación con NaN para no conectar segmentos independientes.
        ksni_xy = np.concatenate([
            ksni_segments[ksni_sel],
            np.full((ksni_sel.sum(), 1, 2), np.nan)
        ], axis=1).reshape(-1, 2)

        ksni_custom = np.repeat(
            ksni_meta[ksni_sel], 3, axis=0
        )

        ksni_fig.add_trace(go.Scattergl(
            x=ksni_xy[:, 0] / 1000,
            y=ksni_xy[:, 1] / 1000,
            mode="lines",
            line=dict(
                color=ksni_color,
                width=ANCHO_CAUCES_KSN
            ),
            customdata=ksni_custom,
            showlegend=False,
            connectgaps=False,
            hovertemplate=(
                "%{customdata[0]}<br>"
                f"ksn: %{{customdata[1]:.2f}} {ksni_unidades}<br>"
                "Distancia desde outlet: "
                "%{customdata[2]:.0f}–%{customdata[3]:.0f} m<br>"
                f"Ventana de {VENTANA_MAPA_KSN_M:g} m: %{{customdata[4]}}"
                "<extra></extra>"
            )
        ))

    # ---------- Fallas como referencias ----------
    for ksni_f in gf_trazas_mapa.itertuples():
        ksni_parts = (
            list(ksni_f.geometry.geoms)
            if hasattr(ksni_f.geometry, "geoms")
            else [ksni_f.geometry]
        )

        for ksni_part in ksni_parts:
            ksni_xy = np.asarray(ksni_part.coords)[:, :2]

            ksni_fig.add_trace(go.Scattergl(
                x=ksni_xy[:, 0] / 1000,
                y=ksni_xy[:, 1] / 1000,
                mode="lines",
                line=dict(
                    color="#5d204a", width=2, dash="dash"
                ),
                name=f"{ksni_f.falla} | {ksni_f.id_traza}",
                hovertemplate=(
                    f"{ksni_f.id_traza}<br>"
                    f"{ksni_f.falla}<extra></extra>"
                )
            ))

    for ksni_f in gf_maule_mapa.itertuples():
        ksni_parts = (
            list(ksni_f.geometry.geoms)
            if hasattr(ksni_f.geometry, "geoms")
            else [ksni_f.geometry]
        )

        for ksni_part in ksni_parts:
            ksni_xy = np.asarray(ksni_part.coords)[:, :2]

            ksni_fig.add_trace(go.Scattergl(
                x=ksni_xy[:, 0] / 1000,
                y=ksni_xy[:, 1] / 1000,
                mode="lines",
                line=dict(
                    color="#555555", width=2, dash="dot"
                ),
                name="Maule: contexto",
                hoverinfo="name"
            ))

    ksni_kp = gf_manuales.loc[
        gf_manuales.ID.isin(RIOS_INTERACTIVOS_KSN)
    ]

    ksni_fig.add_trace(go.Scattergl(
        x=ksni_kp.x_utm / 1000,
        y=ksni_kp.y_utm / 1000,
        mode="markers+text",
        text=ksni_kp.ID,
        textposition="top right",
        marker=dict(
            symbol="star",
            size=11,
            color="white",
            line=dict(color="black", width=1)
        ),
        name="KP manuales",
        hovertext=ksni_kp.id_kp,
        hoverinfo="text"
    ))

    # ---------- Barra: valores reales de ksn ----------
    ksni_ticks = np.unique(np.r_[
        ksni_cmin,
        [
            v for v in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
            if ksni_cmin < v < ksni_cmax
        ],
        ksni_cmax
    ])

    ksni_labels = [f"{v:.1f}" for v in ksni_ticks]
    ksni_labels[0] = "≤" + ksni_labels[0]
    ksni_labels[-1] = "≥" + ksni_labels[-1]

    ksni_point = ksni_segments[
        np.flatnonzero(ksni_ok)[0], 0
    ] / 1000

    ksni_fig.add_trace(go.Scattergl(
        x=[ksni_point[0]] * 2,
        y=[ksni_point[1]] * 2,
        mode="markers",
        marker=dict(
            size=1,
            opacity=0,
            color=[ksni_lmin, ksni_lmax],
            colorscale="Turbo",
            cmin=ksni_lmin,
            cmax=ksni_lmax,
            showscale=True,
            colorbar=dict(
                title=f"ksn ({ksni_unidades})<br>color logarítmico",
                tickvals=np.log10(ksni_ticks),
                ticktext=ksni_labels,
                len=0.7
            )
        ),
        showlegend=False,
        hoverinfo="skip"
    ))

    # ---------- Selector de zoom ----------
    def ksni_ranges(ids):
        xy = np.vstack([
            gf_ksn_perfiles[i][["x_utm", "y_utm"]].to_numpy(float)
            for i in ids
        ]) / 1000

        pad = max(0.5, 0.05 * np.ptp(xy, axis=0).max())

        return {
            "xaxis.range": [
                xy[:, 0].min() - pad,
                xy[:, 0].max() + pad
            ],
            "yaxis.range": [
                xy[:, 1].min() - pad,
                xy[:, 1].max() + pad
            ]
        }

    ksni_views = {
        "Vista general": RIOS_INTERACTIVOS_KSN
    }

    for ksni_label, ksni_ids in {
        "Carranza–Pelluhue": [str(i) for i in range(1, 11)],
        "Sector sur": [str(i) for i in range(8, 19)]
    }.items():
        ksni_ids = [
            i for i in ksni_ids
            if i in RIOS_INTERACTIVOS_KSN
        ]
        if ksni_ids:
            ksni_views[ksni_label] = ksni_ids

    ksni_views.update({
        f"{i}: {gf_nombres[i]}": [i]
        for i in RIOS_INTERACTIVOS_KSN
    })

    ksni_initial = ksni_ranges(RIOS_INTERACTIVOS_KSN)

    ksni_fig.update_layout(
        title=f"ksn local en planta | ventana {VENTANA_MAPA_KSN_M:g} m | theta {THETA_MAPA_KSN:g}",
        height=900,
        dragmode="pan",
        hovermode="closest",
        xaxis=dict(
            title="Este UTM (km)",
            range=ksni_initial["xaxis.range"],
            constrain="domain"
        ),
        yaxis=dict(
            title="Norte UTM (km)",
            range=ksni_initial["yaxis.range"],
            scaleanchor="x",
            scaleratio=1,
            constrain="domain"
        ),
        images=[dict(
            source=ksni_fondo,
            xref="x", yref="y",
            x=ksni_left / 1000,
            y=ksni_top / 1000,
            sizex=(ksni_right - ksni_left) / 1000,
            sizey=(ksni_top - ksni_bottom) / 1000,
            xanchor="left", yanchor="top",
            sizing="stretch",
            opacity=0.5,
            layer="below"
        )],
        updatemenus=[dict(
            buttons=[
                dict(
                    label=label,
                    method="relayout",
                    args=[ksni_ranges(ids)]
                )
                for label, ids in ksni_views.items()
            ],
            x=0, y=1.08,
            xanchor="left", yanchor="top"
        )],
        legend=dict(font=dict(size=10)),
        margin=dict(l=60, r=60, b=50, t=120)
    )


    return ksni_fig
