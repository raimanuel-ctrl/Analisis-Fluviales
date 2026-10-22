
import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import io,base64,html
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LightSource,LogNorm,is_color_like
from matplotlib.lines import Line2D
import rasterio
from rasterio.windows import from_bounds
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
from plotly.offline import get_plotlyjs
import red_fallas as rf

FALLA_COLORS={'Carranza':'#711b53','Pelluhue':'#ac4912','El Yolki':'#49318f'}

def hillshade(dem,bounds,maxside=1400):
    with rasterio.open(dem) as src:
        win=from_bounds(*bounds,src.transform).round_offsets().round_lengths()
        factor=max(1,max(win.width,win.height)/maxside)
        h=max(2,int(win.height/factor));w=max(2,int(win.width/factor))
        z=src.read(1,window=win,out_shape=(h,w),masked=True,boundless=True).astype(float).filled(np.nan)
        valid=np.isfinite(z);fill=np.where(valid,z,np.nanmedian(z))
        shade=LightSource(azdeg=315,altdeg=45).hillshade(fill,dx=abs(src.res[0])*factor,dy=abs(src.res[1])*factor)
        rgba=np.zeros(z.shape+(4,),dtype=np.uint8);rgba[:,:,:3]=(shade[:,:,None]*255).astype('uint8');rgba[:,:,3]=valid*255
        return rgba

def segmentos_xy(edges):
    return np.array([np.array(g.coords)[[0,-1],:2]/1000 for g in edges.geometry if g.geom_type=='LineString'])

def plot_mapa(edges,fallas,poly,crossings,dem,out,window=1000.):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    edges=edges[edges.geom_type.eq('LineString')].copy()
    good=edges.ksn.gt(0)&np.isfinite(edges.ksn)
    lo,hi=np.nanpercentile(edges.loc[good,'ksn'],[5,95]);lo=max(lo,1e-3);hi=max(hi,lo*1.01)
    norm=LogNorm(lo,hi,clip=True);segs=segmentos_xy(edges)
    systems=fallas.nombre_base.unique().tolist()
    fig,axes=plt.subplots(1,1+min(3,len(systems)),figsize=(6+5*min(3,len(systems)),8),layout='constrained')
    axes=np.atleast_1d(axes);boxes=[poly.bounds]
    for name in systems[:3]:
        b=fallas[fallas.nombre_base.eq(name)].total_bounds;pad=2500
        boxes.append((b[0]-pad,b[1]-pad,b[2]+pad,b[3]+pad))
    collection=None
    for i,(ax,bounds) in enumerate(zip(axes,boxes)):
        bg=hillshade(dem,bounds,maxside=1000)
        ax.imshow(bg,extent=np.array(bounds)[[0,2,1,3]]/1000,origin='upper',alpha=.6,zorder=0)
        ax.add_collection(LineCollection(segs[~good.to_numpy()],colors='#8f959c',linewidths=.6,zorder=1))
        collection=LineCollection(segs[good.to_numpy()],array=edges.loc[good,'ksn'].to_numpy(),cmap='turbo',norm=norm,linewidths=.75,zorder=2)
        ax.add_collection(collection)
        for _,f in fallas.iterrows():
            lines=list(f.geometry.geoms) if f.geometry.geom_type=='MultiLineString' else [f.geometry]
            for line in lines:
                xy=np.array(line.coords)/1000;ax.plot(xy[:,0],xy[:,1],color='white',lw=3.5,zorder=3)
                ax.plot(xy[:,0],xy[:,1],color=FALLA_COLORS.get(f.nombre_base,'black'),lw=1.6,ls='--',zorder=4)
        selected=crossings[crossings.get('estado',pd.Series(index=crossings.index,dtype=str)).eq('seleccionado')]
        if len(selected):
            ax.scatter(selected.x_cruce/1000,selected.y_cruce/1000,c='#101010',edgecolors='white',s=25,zorder=5)
            for _,v in selected.drop_duplicates('perfil_id').iterrows():ax.annotate(v.perfil_id,(v.x_cruce/1000,v.y_cruce/1000),xytext=(3,4),textcoords='offset points',fontsize=7,zorder=6)
        ax.set_xlim(bounds[0]/1000,bounds[2]/1000);ax.set_ylim(bounds[1]/1000,bounds[3]/1000)
        ax.set_aspect('equal');ax.set_xlabel('Este UTM (km)');ax.set_ylabel('Norte UTM (km)');ax.grid(alpha=.15)
        ax.set_title('Red completa del polígono' if i==0 else systems[i-1],fontsize=12)
    fig.colorbar(collection,ax=axes.tolist(),label=r'$k_{sn}$ (m$^{0.8}$), colores logarítmicos · percentiles 5–95',shrink=.65,pad=.02)
    fig.suptitle(f'ksn local · ventana {window:g} m · θ=0,40 · todos los cauces con A ≥ {edges.area_km2.min():.2f} km²\nNegro: cruces con perfiles seleccionados; gris: ksn no válido',fontsize=13)
    target=out/'Mapa_red_ksn_fallas.png';fig.savefig(target,dpi=160,facecolor='white');plt.close(fig)
    # HTML con toda la geometría de la red; hover submuestreado para mantener fluidez.
    f=go.Figure();val=edges.ksn.to_numpy();finite=np.isfinite(val)&(val>0)
    classes=np.clip(np.floor((np.log10(np.maximum(np.nan_to_num(val,nan=lo),lo))-np.log10(lo))/(np.log10(hi)-np.log10(lo))*24),0,23).astype(int)
    for k in range(24):
        coords=segs[finite&(classes==k)]
        if not len(coords):continue
        xx=np.c_[coords[:,:,0],np.full(len(coords),np.nan)].ravel();yy=np.c_[coords[:,:,1],np.full(len(coords),np.nan)].ravel()
        color=sample_colorscale('Turbo',[k/23])[0]
        f.add_trace(go.Scattergl(x=np.round(xx,3).tolist(),y=np.round(yy,3).tolist(),mode='lines',line={'color':color,'width':1.6},hoverinfo='skip',showlegend=False,name='ksn'))
    invalid=segs[~finite]
    if len(invalid):
        f.add_trace(go.Scattergl(x=np.c_[invalid[:,:,0],np.full(len(invalid),np.nan)].ravel(),
            y=np.c_[invalid[:,:,1],np.full(len(invalid),np.nan)].ravel(),mode='lines',line={'color':'#999999','width':1},name='ksn no válido',hoverinfo='skip'))
    sample=np.flatnonzero(finite)[::8];center=segs[sample].mean(axis=1);q=edges.iloc[sample]
    f.add_trace(go.Scattergl(x=center[:,0],y=center[:,1],mode='markers',showlegend=False,
        marker={'size':3,'opacity':.1,'color':np.log10(q.ksn),'cmin':np.log10(lo),'cmax':np.log10(hi),'colorscale':'Turbo','showscale':False},
        customdata=q[['ksn','R2','area_km2','orden','span','gid_u']].to_numpy(),
        hovertemplate='ksn %{customdata[0]:.2f} m^0.8<br>R² %{customdata[1]:.3f}<br>Área %{customdata[2]:.2f} km²<br>Orden %{customdata[3]:.0f}<br>Ventana efectiva %{customdata[4]:.0f} m<br>Nodo %{customdata[5]:.0f}<extra></extra>'))
    for name in systems:
        xx=[];yy=[]
        for geom in fallas[fallas.nombre_base.eq(name)].geometry:
            for line in (list(geom.geoms) if geom.geom_type=='MultiLineString' else [geom]):
                a=np.array(line.coords)/1000;xx.extend(a[:,0].tolist()+[None]);yy.extend(a[:,1].tolist()+[None])
        f.add_trace(go.Scatter(x=xx,y=yy,mode='lines',line={'color':'white','width':5},hoverinfo='skip',showlegend=False))
        f.add_trace(go.Scatter(x=xx,y=yy,mode='lines',name=name,line={'color':FALLA_COLORS.get(name,'black'),'width':2.2,'dash':'dash'}))
    if len(selected):f.add_trace(go.Scatter(x=selected.x_cruce/1000,y=selected.y_cruce/1000,mode='markers+text',text=selected.perfil_id,textposition='top center',marker={'color':'black','size':8,'line':{'color':'white','width':1}},name='Perfiles seleccionados',hovertemplate='%{text}<extra>Cruce cartográfico</extra>'))
    ticks=np.geomspace(lo,hi,6)
    f.add_trace(go.Scatter(x=[None],y=[None],mode='markers',showlegend=False,
        marker={'color':[np.log10(lo)],'cmin':np.log10(lo),'cmax':np.log10(hi),'colorscale':'Turbo','showscale':True,
                'colorbar':{'title':'ksn (m^0.8)','tickvals':np.log10(ticks),'ticktext':[f'{x:.1f}' for x in ticks]}}))
    bounds=poly.bounds;bg=hillshade(dem,bounds)
    from PIL import Image
    buf=io.BytesIO();Image.fromarray(bg).save(buf,format='PNG')
    image='data:image/png;base64,'+base64.b64encode(buf.getvalue()).decode()
    f.add_layout_image({'source':image,'xref':'x','yref':'y','x':bounds[0]/1000,'y':bounds[3]/1000,
       'sizex':(bounds[2]-bounds[0])/1000,'sizey':(bounds[3]-bounds[1])/1000,'layer':'below','sizing':'stretch','opacity':.65})
    f.update_layout(template='plotly_white',height=1000,title=f'Red completa y ksn · ventana {window:g} m · zoom con rueda',
        xaxis={'title':'Este UTM (km)','range':[bounds[0]/1000,bounds[2]/1000]},
        yaxis={'title':'Norte UTM (km)','range':[bounds[1]/1000,bounds[3]/1000],'scaleanchor':'x'},legend={'orientation':'h','y':-.1})
    htmlpath=out/'Mapa_red_ksn_interactivo.html';f.write_html(htmlpath,include_plotlyjs=True,config={'scrollZoom':True,'displaylogo':False})
    return target,htmlpath,{'ksn_color_min':float(lo),'ksn_color_max':float(hi),'n_segmentos':len(edges),'hover_cada_n_segmentos':8}

def plot_perfiles(profiles,crossings,out,window=1000.):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);paths=[];parts=[];summary=[]
    for pid,p in profiles.items():
        stats=rf.ventana_ols(p.distancia_m,p.chi_m,p.elevacion_m,W=window)
        p['ksn_perfil']=stats['ksn'];p['R2_ksn']=stats['R2']
        codes=p.codigo_lito.drop_duplicates().tolist()
        colors={c:p[p.codigo_lito.eq(c)].color_lito.iloc[0] for c in codes}
        colors={k:v if is_color_like(v) else '#999999' for k,v in colors.items()}
        fig,(ax,ak)=plt.subplots(2,1,figsize=(13,8),sharex=True,gridspec_kw={'height_ratios':[2.3,1]},layout='constrained')
        points=np.c_[p.chi_m,p.elevacion_m];segments=np.stack([points[:-1],points[1:]],axis=1)
        ax.add_collection(LineCollection(segments,colors=[colors[v] for v in p.codigo_lito.iloc[:-1]],linewidths=2.5))
        ax.autoscale();ak.plot(p.chi_m,p.ksn_perfil,color='#20587b',lw=1)
        cr=crossings[crossings.perfil_id.eq(pid)] if len(crossings) else crossings
        for _,v in cr.iterrows():
            for axis in [ax,ak]:axis.axvline(v.chi_cruce_m,color=FALLA_COLORS.get(v.nombre_base,'black'),ls='--',lw=1)
        for _,v in cr.drop_duplicates('nombre_base').iterrows():
            ax.annotate(v.nombre_base,(v.chi_cruce_m,1),xycoords=('data','axes fraction'),xytext=(3,-3),textcoords='offset points',rotation=90,va='top',fontsize=8,color=FALLA_COLORS.get(v.nombre_base,'black'))
        ax.set_ylabel('Elevación DEM (m)');ak.set_ylabel('ksn (m^0.8)');ak.set_xlabel('χ (m), desde salida D8 de este cauce')
        ax.grid(alpha=.18);ak.grid(alpha=.18);ax.set_xlim(p.chi_m.min(),p.chi_m.max());ak.set_yscale('log')
        legends=[]
        for c in codes:
            desc=p[p.codigo_lito.eq(c)].descripcion_lito.iloc[0]
            label=c if c in ['Sin cobertura','Superposición'] else c+' · '+str(desc)[:72]
            legends.append(Line2D([],[],color=colors[c],lw=3,label=label))
        fig.legend(handles=legends,loc='outside lower center',ncol=2,fontsize=8)
        missing=float(p.codigo_lito.eq('Sin cobertura').mean())
        fig.suptitle(f'{pid} · perfil χ con litología y cruces de falla\nLongitud {p.distancia_m.iloc[-1]/1000:.1f} km · sin cobertura litológica: {missing:.0%} · ksn ventana {window:g} m',fontsize=13)
        target=out/f'{pid}_chi_litologia.png';fig.savefig(target,dpi=165,facecolor='white');plt.close(fig);paths.append(target)
        # Perfil interactivo: colores de unidad, χ exacto de cruces y hover con fuente.
        f=go.Figure()
        for code in codes:
            xx=[];yy=[]
            for j in range(len(p)-1):
                if p.codigo_lito.iloc[j]==code:xx.extend([float(p.chi_m.iloc[j]),float(p.chi_m.iloc[j+1]),None]);yy.extend([float(p.elevacion_m.iloc[j]),float(p.elevacion_m.iloc[j+1]),None])
            f.add_trace(go.Scatter(x=xx,y=yy,mode='lines',line={'color':colors[code],'width':3},name=code,hoverinfo='skip'))
        f.add_trace(go.Scatter(x=p.chi_m.tolist(),y=p.elevacion_m.tolist(),mode='markers',showlegend=False,
            marker={'size':4,'opacity':.05,'color':'#222'},customdata=p[['codigo_lito','descripcion_lito','fuentes_lito','area_m2','ksn_perfil']].to_numpy().tolist(),
            hovertemplate='χ %{x:.2f} m<br>Z %{y:.2f} m<br>%{customdata[0]} · %{customdata[1]}<br>Fuente: %{customdata[2]}<br>Área %{customdata[3]:.0f} m²<br>ksn %{customdata[4]:.2f}<extra></extra>'))
        for _,v in cr.iterrows():f.add_vline(x=v.chi_cruce_m,line_color=FALLA_COLORS.get(v.nombre_base,'black'),line_dash='dash')
        f.update_layout(template='plotly_white',height=540,title=pid+' · χ, litología y fallas',xaxis_title='χ (m)',yaxis_title='Elevación DEM (m)',legend={'orientation':'h','y':-.22},margin={'b':150})
        block=f.to_html(include_plotlyjs=False,full_html=False,div_id='perfil-'+pid,config={'displaylogo':False,'scrollZoom':True})
        parts.append('<section id="'+pid+'"><h2>'+pid+'</h2><p>Sin cobertura litológica: '+f'{missing:.0%}'+'. Gris significa ausencia de cartografía, no uniformidad de roca. Contactos cartográficos indicativos.</p>'+block+f'<a href="{target.name}">Figura PNG con ksn</a></section>')
        p.to_csv(out/f'{pid}_perfil_completo.csv',index=False,encoding='utf-8-sig')
        summary.append({'perfil_id':pid,'longitud_m':float(p.distancia_m.iloc[-1]),'n_nodos':len(p),'n_cruces':len(cr),
                        'fraccion_sin_litologia':missing,'unidades':' / '.join(codes)})
    doc='<!doctype html><html lang="es"><meta charset="utf-8"><title>Perfiles χ, litología y fallas</title><style>body{font:15px system-ui;color:#183247;max-width:1450px;margin:25px auto;padding:0 20px}section{margin:40px 0}nav{display:flex;gap:15px;flex-wrap:wrap}</style><h1>Perfiles que cruzan fallas · litología cartográfica</h1><p>Perfiles completos dentro de la red umbralizada, sin cortarlos al borde del polígono. Los tramos compartidos no son observaciones independientes. Una coincidencia espacial no demuestra actividad por sí sola.</p><nav>'+' '.join(f'<a href="#{p}">{p}</a>' for p in profiles)+'</nav><script>'+get_plotlyjs()+'</script>'+''.join(parts)+'</html>'
    htmlpath=out/'Perfiles_chi_litologia_interactivos.html';htmlpath.write_text(doc,encoding='utf8')
    return paths,htmlpath,pd.DataFrame(summary)
