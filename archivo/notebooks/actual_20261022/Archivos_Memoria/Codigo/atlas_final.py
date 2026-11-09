"""Atlas a partir de tablas calculadas; no ejecuta nada al importar."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from functools import lru_cache
from pathlib import Path
import textwrap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


@lru_cache(maxsize=4)
def capa(path):
    import geopandas as gpd
    g=gpd.read_file(path)
    if g.crs is None:raise ValueError('Capa sin CRS: '+str(path))
    return g.to_crs('EPSG:32718')


def number(x,scientific=False):
    try:
        x=float(x)
        return ('—' if not np.isfinite(x) else f'{x:.3e}' if scientific else f'{x:.3f}')
    except (TypeError,ValueError):return '—'


def atlas_rio(rid,resultados,perfiles,cruces,segmentos,asociaciones,cfg):
    r=resultados.loc[resultados.rid.eq(rid)].iloc[0]
    p=perfiles[rid];df=p['df'];s,z,chi=p['s'],p['z'],p['chi']
    q=segmentos.loc[segmentos.ID.astype(str).eq(str(rid))]
    t=asociaciones.loc[asociaciones.ID.astype(str).eq(str(rid))]
    fig=plt.figure(figsize=(23,16),layout='constrained')
    gs=fig.add_gridspec(4,3,width_ratios=[1.12,1.1,1.4],height_ratios=[1,1,.8,.8])
    long=fig.add_subplot(gs[0,0]);zoom=fig.add_subplot(gs[0,1]);ch=fig.add_subplot(gs[1,0])
    res=fig.add_subplot(gs[1,1]);ks=fig.add_subplot(gs[2,:2]);notes=fig.add_subplot(gs[3,:2])
    mp=fig.add_subplot(gs[:2,2]);tb=fig.add_subplot(gs[2:,2])
    long.plot(s/1000,z,color='#16465c',lw=1)
    zoom.plot(s/1000,z,color='#16465c',lw=1)
    ch.plot(chi/1000,z,color='#16465c',lw=.9)
    if len(q):
        for ax,x0,x1 in [(long,q.s0/1000,q.s1/1000),(ch,q.chi0/1000,q.chi1/1000)]:
            lines=np.stack([np.c_[x0,q.z0],np.c_[x1,q.z1]],axis=1)
            cols=[v if isinstance(v,str) and v.startswith('#') else '#aaa' for v in q.color_hex]
            ax.add_collection(LineCollection(lines,colors=cols,linewidths=2))
    c=cruces.loc[cruces.rid.eq(rid)]
    if np.isfinite(r.H_m):
        for ax in [long,zoom,ch]:ax.axhline(r.H_m,ls='--',lw=.8,color='#af702e')
        for _,cr in c.iterrows():
            long.scatter(cr.s_PT_m/1000,r.H_m,s=20,facecolor='white',edgecolor='#af702e')
            zoom.scatter(cr.s_PT_m/1000,r.H_m,s=20,facecolor='white',edgecolor='#af702e')
            ch.scatter(cr.chi_PT_m/1000,r.H_m,s=20,facecolor='white',edgecolor='#af702e')
    if np.isfinite(r.s_PT_m):
        for ax,x in [(long,r.s_PT_m/1000),(zoom,r.s_PT_m/1000),(ch,r.chi_PT_m/1000)]:
            ax.scatter(x,r.H_m,s=42,color='#af702e',zorder=7)
        end=min(s[-1],max(500,r.s_PT_m*1.25,float(c.s_PT_m.max())*1.1 if len(c) else 0))
        sel=s<=end
        zz=z[sel];bottom=min(np.min(zz),r.H_m);top=max(np.max(zz),r.H_m)
        zoom.set_xlim(0,end/1000);zoom.set_ylim(bottom-max(3,(top-bottom)*.1),top+max(3,(top-bottom)*.1))
        mask=s<=r.s_PT_m
        long.axvspan(0,r.s_PT_m/1000,color='#4e87a0',alpha=.12)
        if np.isfinite(r.Mchi):
            fitted=r.intercepto_m+r.Mchi*chi[mask]
            ch.plot(chi[mask]/1000,fitted,color='#b64533',lw=1.4,label='OLS del tramo')
            res.plot(chi[mask]/1000,z[mask]-fitted,lw=.9,color='#16465c')
    else:
        zoom.text(.5,.5,str(r.estado_K),ha='center',va='center',transform=zoom.transAxes,wrap=True)
    if np.isfinite(r.s_KP_m):
        k=np.interp(r.s_KP_m,s,z)
        long.scatter(r.s_KP_m/1000,k,marker='*',color='black',s=90,zorder=8,label='KP manual')
        ch.scatter(r.chi_KP_m/1000,k,marker='*',color='black',s=90,zorder=8)
    ks.plot(s/1000,df.ksn,lw=.9,color='#16465c')
    if np.isfinite(r.s_KP_m):ks.axvline(r.s_KP_m/1000,lw=1,color='black',ls=':')
    res.axhline(0,color='#888',lw=.7)
    long.set(title='a · Longitudinal y litología',xlabel='Distancia desde outlet (km)',ylabel='Elevación (m)')
    zoom.set(title='b · Tramo inferior y cruces de cota',xlabel='Distancia desde outlet (km)',ylabel='Elevación (m)')
    ch.set(title=f'c · χ–z, θ={cfg["theta"]:g}',xlabel='χ (km)',ylabel='Elevación (m)')
    res.set(title=f'd · Residuos · R²={number(r.R2)} · RMSE={number(r.RMSE_m)} m',xlabel='χ del tramo (km)',ylabel='Residuo (m)')
    ks.set(title=f'e · ksn local, ventana {cfg["ventana_m"]:g} m',xlabel='Distancia desde outlet (km)',ylabel=f'ksn (m^{2*cfg["theta"]:g})')
    for ax in [long,zoom,ch,res,ks]:ax.grid(alpha=.14)

    g=capa(cfg['geologia_integrada']);f=capa(cfg['chaf'])
    xmin,xmax=df.x_utm.min()-1000,df.x_utm.max()+1000
    ymin,ymax=df.y_utm.min()-1000,df.y_utm.max()+1000
    if len(t):
        xmin=min(xmin,t['Este UTM18S (m)'].min()-500);xmax=max(xmax,t['Este UTM18S (m)'].max()+500)
        ymin=min(ymin,t['Norte UTM18S (m)'].min()-500);ymax=max(ymax,t['Norte UTM18S (m)'].max()+500)
    g=g.cx[xmin:xmax,ymin:ymax];f=f.cx[xmin:xmax,ymin:ymax]
    if len(g):g.plot(ax=mp,color=g.color_hex.fillna('#ddd').tolist(),edgecolor='#888',linewidth=.15,alpha=.8,rasterized=True)
    if len(f):f.plot(ax=mp,color='#a53b40',lw=1,ls='--')
    mp.plot(df.x_utm,df.y_utm,color='#183546',lw=1.4)
    mp.scatter(df.x_utm.iloc[0],df.y_utm.iloc[0],marker='s',color='#183546',s=35,label='Outlet')
    for _,tr in t.iterrows():
        mp.scatter(tr['Este UTM18S (m)'],tr['Norte UTM18S (m)'],marker='D',s=35,color='#af702e',zorder=7)
        mp.plot([tr['Este UTM18S (m)'],tr.x_proyectado],[tr['Norte UTM18S (m)'],tr.y_proyectado],ls=':',lw=.6,color='#af702e')
    if np.isfinite(r.s_PT_m):
        mp.scatter(np.interp(r.s_PT_m,s,df.x_utm),np.interp(r.s_PT_m,s,df.y_utm),s=45,color='#af702e',edgecolor='white',zorder=8)
    if np.isfinite(r.s_KP_m):
        mp.scatter(np.interp(r.s_KP_m,s,df.x_utm),np.interp(r.s_KP_m,s,df.y_utm),marker='*',s=110,color='black',zorder=9)
    mp.set(xlim=(xmin,xmax),ylim=(ymin,ymax),title='f · Contexto en planta · UTM 18S',xlabel='Este (m)',ylabel='Norte (m)')
    mp.set_aspect('equal');mp.ticklabel_format(style='plain',useOffset=False)
    span=xmax-xmin;l=span*.2;x=xmin+span*.08;y=ymin+(ymax-ymin)*.06
    mp.plot([x,x+l],[y,y],color='black',lw=3);mp.text(x+l/2,y,f'{l/1000:.1f} km',ha='center',va='bottom',fontsize=8)
    mp.annotate('N',xy=(.92,.96),xytext=(.92,.83),xycoords='axes fraction',arrowprops={'arrowstyle':'->'},ha='center')
    labels=g[['grupo_lito','color_hex']].drop_duplicates() if len(g) else None
    handles=[] if labels is None else [Patch(facecolor=v.color_hex if isinstance(v.color_hex,str) else '#ddd',label=str(v.grupo_lito)) for _,v in labels.iterrows()]
    handles.extend([Line2D([0],[0],marker='D',ls='',color='#af702e',label='Terraza XY'),
                    Line2D([0],[0],marker='*',ls='',color='black',label='KP manual')])
    mp.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.01,1),fontsize=7,frameon=False)

    tb.axis('off')
    rows=[['K viaje χ',number(r.K_George,True),number(r.error_K_George,True)],
          ['K pendiente',number(r.K_pendiente,True),number(r.error_K_pendiente,True)],
          ['K con U publicada',number(r.K_JM_pendiente,True),number(r.error_K_JM_pendiente,True)],
          ['K regional transferido',number(r.get('K_regional'),True),number(r.get('error_K_regional'),True)],
          ['U terraza (mm/a)',number(r.U_terraza),number(r.error_U_terraza)],
          ['U fluvial condicionado',number(r.U_George),number(r.error_U_George)],
          ['U publicada (mm/a)',number(r.U_JM_m_ka),number(r.error_U_JM_m_ka)],
          ['T (ka)',number(r.get('edad_ka')),number(r.get('error_edad_ka'))],
          ['H (m)',number(r.H_m),number(r.error_H_publicado_m)],
          ['Cruces / elegido',str(int(r.n_cruces)),number(r.cruce_elegido)],
          ['Tiempo viaje KP (ka)',number(r.t_viaje_KP_ka),'condicional']]
    table=tb.table(cellText=rows,colLabels=['g · Magnitud','Valor','Error publicado\npropagado'],
                   loc='center',cellLoc='left',colWidths=[.43,.24,.33])
    table.auto_set_font_size(False);table.set_fontsize(8);table.scale(1,1.8)
    for (i,j),ce in table.get_celld().items():
        ce.set_edgecolor('#d9e1e6')
        if i==0:ce.set_facecolor('#224c63');ce.set_text_props(color='white',weight='bold',fontsize=8)
        elif i%2:ce.set_facecolor('#f2f6f8')
    notes.axis('off')
    text=(f'Terraza: {r.terraza_id or "sin selección"} · MIS: {r.MIS or "sin dato"} · KP: {r.id_KP or "sin selección"}. '
          f'Estado: {r.estado_K}. Estado pendiente: {r.estado_K_pendiente}. '
          f'K tiene unidades m^{1-2*cfg["theta"]:g}/año para n=1. B={r.chi_B_m:g} m de χ; DEM fijo. '
          'Las barras propagan magnitudes publicadas; no constituyen un intervalo de confianza común ni incluyen incertidumbre de asociación, '
          'paleocosta o elección entre cruces. K regional usa U transferido y no sustituye K local. '
          'El tiempo de viaje es condicionado; no es una edad tectónica. La litología proviene de cartografía y conserva vacíos y ambigüedades.')
    notes.text(0,1,'\n'.join(textwrap.wrap(text,155)),va='top',fontsize=9,linespacing=1.5,transform=notes.transAxes)
    fig.suptitle(f'{rid:02d} · {r.rio} · Atlas de calibración y contexto',fontsize=20,weight='bold')
    return fig
