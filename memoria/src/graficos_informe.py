"""Regenerate the thesis graphs. Cartographic compositions are excluded.

The preserved scientific values are replotted; this does not execute DEM routing.
Figure 2 is a bibliographic reproduction. Figures 16 and 29 need local comparison
inputs, which are intentionally absent from the public distribution.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,hashlib,shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from . import analisis as a

DATA=a.DATA/'graficos'
MAPS={1,3,12,17,25,27,32,40}
BLUE,GREEN,RED,GRAY='#185984','#28836b','#c34360','#9aa5a9'
PALETTE=plt.colormaps['tab20'](np.linspace(0,1,18))


class Graficos:
    def __init__(self,output,privadas=None,formats=('png','svg')):
        self.out=Path(output);self.out.mkdir(parents=True,exist_ok=False)
        self.private=Path(privadas) if privadas else a.BASE/'entradas_privadas'
        self.formats=formats;self.k=a.leer('catalogo_kp.csv').set_index('rid')
        self.p={r:a.perfil(r) for r in range(1,19)}
        self.t=a.tiempos_finales();self.u=a.alzamiento_tramo();self.means=a.calibraciones_centrales()
        self.cal=a.leer('calibraciones_k_por_pt.csv');self.orig=a.leer('origenes_temporales.csv').set_index('rid')
        self.lat=self.read('latitudes_kp.csv').set_index('rid').latitud_sur
        self.cuts=a.leer('reloca_fronteras.csv');self.c2=float(self.cuts.query("frontera=='C2'").chi_m.iloc[0])
        self.rates=a.incision_ventanas();self.record=[]

    def read(self,name):return pd.read_csv(DATA/name,float_precision='round_trip')

    def pair(self,rows=1):return plt.subplots(rows,2,figsize=(8,3.7*rows),layout='constrained',squeeze=False)

    def star(self,ax,x,y,r,label=None):
        ax.scatter(x,y,marker='*',s=65,facecolors='white' if r in [11,14,17] else RED,edgecolors=RED,zorder=6,label=label)

    def point(self,ax,x,y,r,annotate=True):
        self.star(ax,x,y,r)
        if annotate:ax.annotate(f'R{r:02d}',(x,y),xytext=(4,5 if r%2 else -12),textcoords='offset points',fontsize=7)

    def curve(self,ax,x,y,**kw):
        x=np.asarray(x);y=np.asarray(y);good=x>=-1e-8
        if (~good).any():ax.plot(x[~good],y[~good],color=GRAY,ls=':',lw=.8)
        ax.plot(x[good],y[good],**kw);ax.axvline(0,color=GRAY,lw=.6,ls='--')

    def timecurve(self,r,mis,method='G'):
        v=self.t.query('rid==@r and mis==@mis').iloc[0]
        return a.tiempo(self.p[r].chi_m,v.chi_origen_m,v[f'K_{method}_inferior_m02_anio'],v[f'K_{method}_superior_m02_anio'],self.c2 if r==6 else None),v

    def save(self,fig,n,notes=''):
        for ax in fig.axes:
            if ax.axison:ax.grid(alpha=.14)
        name=f'Figura_{n:02d}'+('_solo_perfil' if n==14 else '')
        outputs=[]
        for ext in self.formats:
            p=self.out/f'{name}.{ext}';fig.savefig(p,dpi=200,facecolor='white');outputs.append(p.name)
        plt.close(fig);self.record.append(dict(figura=n,estado='generada',archivos=outputs,nota=notes))

    def conceptual(self,n):
        if n==4:
            fig,ax=self.pair();x=np.linspace(0,1,200);z=np.where(x<.55,x**.6*.7,.7*.55**.6+1.1*(x-.55))
            for panel,horizontal,title in [(ax[0,0],x,'Perfil longitudinal'),(ax[0,1],np.where(x<.55,.8*z,.8*z[109]+1.3*(z-z[109])),'Perfil χ–elevación')]:
                panel.plot(horizontal,z,c=BLUE);panel.scatter(horizontal[109],z[109],c=BLUE);panel.set(title=title,xlabel='Coordenada relativa',ylabel='Elevación relativa')
        elif n==5:
            fig,ax=plt.subplots(figsize=(7,4),layout='constrained');x=np.linspace(0,1,200)
            for front,col,style,lab in [(0,GRAY,'--','Perfil previo'),(.55,BLUE,'-','Estado intermedio'),(.8,GREEN,'-','Estado posterior')]:
                inherited=.35*x+.25
                z=np.where(x<front,x*(.35*front+.25)/front,inherited) if front else inherited
                ax.plot(x,z,c=col,ls=style,label=lab)
            ax.annotate('Propagación aguas arriba',xy=(.8,.6),xytext=(.3,.6),arrowprops={'arrowstyle':'->'});ax.set(xlabel='Distancia relativa',ylabel='Elevación relativa');ax.legend()
        elif n==6:
            fig,axes=self.pair()
            for ax,text in zip(axes.flat,[r'$K_G=\Delta\chi/(A_0^\theta T)$',r'$K_P=U_{terraza}/k_{sn}$'+'\n'+r'$k_{sn}=M_\chi A_0^\theta$']):
                ax.plot([0,.5,1],[0,.35,1],c=BLUE);ax.axhline(.35,c=GREEN,ls=':');ax.text(.05,.82,text,transform=ax.transAxes);ax.scatter(0,0,marker='v',c=BLUE);ax.scatter(.5,.35,c=BLUE);ax.annotate('B',(0,0),xytext=(6,5),textcoords='offset points');ax.annotate('Pₜ',(.5,.35),xytext=(5,6),textcoords='offset points');ax.text(.02,.38,'H',transform=ax.transAxes);ax.set(xlabel='χ relativa',ylabel='Elevación relativa')
        else:
            fig,ax=plt.subplots(figsize=(8,6));ax.set_axis_off()
            boxes=[(.04,.80,'DEM, red y perfiles'),(.57,.80,'Terrazas y cronología MIS'),(.04,.57,'Revisión de KP y ajustes'),(.57,.57,'B–Pₜ y calibración de K'),(.30,.29,'KP adoptados, K y recorridos'),(.30,.06,'Tiempos, tasas y comparación regional')]
            for x,y,text in boxes:ax.add_patch(FancyBboxPatch((x,y),.37,.12,boxstyle='round,pad=.015',facecolor='white',edgecolor=BLUE));ax.text(x+.185,y+.06,text,ha='center',va='center',fontsize=10)
            for xy,xytext in [((.225,.70),(.225,.80)),((.755,.70),(.755,.80)),((.48,.42),(.225,.57)),((.48,.42),(.755,.57)),((.48,.20),(.48,.29))]:ax.annotate('',xy=xy,xytext=xytext,arrowprops={'arrowstyle':'->','color':BLUE})
        self.save(fig,n,'Esquema conceptual, sin datos medidos; geometría ilustrativa.')

    def figura_8(self):
        fig,axes=self.pair(2)
        for ax,r in zip(axes.flat,[1,4,11,14]):
            c=self.read(f'R{r:02d}_perfil_costa.csv');p=self.p[r];k=self.k.loc[r]
            delta=np.hypot(c.x_utm-p.x_utm.iloc[0],c.y_utm-p.y_utm.iloc[0]);idx=delta.idxmin();assert delta.loc[idx]<.001
            start=c.loc[idx,'distancia_costa_m'];sk=k.s_m-p.distancia_m.iloc[0]+start
            ax.plot(c.distancia_costa_m/1000,c.elevacion_m,c=BLUE);ax.axvspan(0,start/1000,color='#e7edf2');ax.axvline(start/1000,c=BLUE,ls='--');ax.scatter(start/1000,p.elevacion_m.iloc[0],marker='s',c=BLUE)
            self.star(ax,sk/1000,k.z_m,r);ax.set(title=f'R{r:02d}',xlabel='Distancia desde costa (km)',ylabel='Elevación DEM (m)')
        self.save(fig,8)

    def figura_9(self):
        fig,ax=plt.subplots(figsize=(8,4),layout='constrained');p=self.p[9];c=self.orig.loc[9]
        ax.plot(p.chi_m/1000,p.elevacion_m,c=BLUE);ax.axvspan(0,c.chi_origen_m/1000,color='#e7edf2');ax.axvline(c.chi_origen_m/1000,c=GRAY,ls='--');self.star(ax,self.k.loc[9].chi_m/1000,self.k.loc[9].z_m,9);ax.set(xlabel='χ desde outlet (km)',ylabel='Elevación DEM (m)');self.save(fig,9)

    def figura_10(self):
        fig,ax=plt.subplots(figsize=(7,5),layout='constrained');p=self.p[5];q=self.cal.query('rid==5');lo=q.s_b_m.min()-300;hi=q.s_pt_m.max()+500;s=p.distancia_m.between(lo,hi)
        ax.plot(p.distancia_m[s]/1000,p.elevacion_m[s],c=BLUE);b=q.iloc[0];ax.scatter(b.s_b_m/1000,b.z_b_m,marker='v',c=BLUE,label='B')
        for j,v in enumerate(q.itertuples()):ax.scatter(v.s_pt_m/1000,v.H_dem_m,label=v.etiqueta_pt,marker=['o','s','D'][j],s=35)
        ax.axhline(b.H_dem_m,c=GRAY,ls=':');ax.set(xlabel='Distancia del perfil (km)',ylabel='Elevación DEM (m)');ax.legend();self.save(fig,10)

    def figura_11(self):
        fig,axes=plt.subplots(2,1,figsize=(9,7),layout='constrained',sharey=True)
        for ax,rs in zip(axes,[range(1,10),range(10,19)]):
            for j,r in enumerate(rs):
                p=self.p[r];ax.plot((p.chi_m-p.chi_m.iloc[0]+j*650)/1000,p.elevacion_m-p.elevacion_m.iloc[0],c=PALETTE[r-1],label=f'R{r:02d}')
                if r in self.k.index:self.star(ax,(self.k.loc[r].chi_m-p.chi_m.iloc[0]+j*650)/1000,self.k.loc[r].z_m-p.elevacion_m.iloc[0],r)
            ax.set(xlabel='χ trasladado horizontalmente (km)',ylabel='Elevación sobre outlet (m)');ax.legend(ncol=5,fontsize=7)
        self.save(fig,11,'Desplazamiento gráfico de 650 m entre curvas; sin significado geográfico.')

    def figura_13(self):
        fig,axes=self.pair(2)
        for row,r in enumerate([7,9]):
            q=self.cal.query('rid==@r and etiqueta_pt=="intermedio_real"');p=self.p[r];mask=p.distancia_m<=q.s_pt_m.max()*1.08
            for ax,col in zip(axes[row],['distancia_m','chi_m']):
                ax.plot(p[col][mask]/1000,p.elevacion_m[mask],c=BLUE)
                for v in q.itertuples():
                    b=v.s_b_m if col=='distancia_m' else v.chi_b_m;pt=v.s_pt_m if col=='distancia_m' else v.chi_pt_m;colour=BLUE if v.mis=='5c' else GREEN
                    ax.axhline(v.H_dem_m,c=colour,ls=':');ax.scatter(b/1000,v.z_b_m,c=colour,marker='v');ax.scatter(pt/1000,v.H_dem_m,c=colour,label='MIS '+v.mis)
                    if col=='chi_m':ax.annotate(f'Δχ={pt-b:.0f} m',((pt+b)/2000,v.H_dem_m-15),fontsize=8,ha='center')
                ax.set(title=f'R{r:02d}',xlabel='Distancia (km)' if col=='distancia_m' else 'χ (km)',ylabel='Elevación DEM (m)');ax.legend(fontsize=8)
        self.save(fig,13)

    def figura_14(self):
        fig,ax=plt.subplots(figsize=(9,3),layout='constrained');p=self.p[9];ax.plot(p.distancia_m/1000,p.elevacion_m,c=BLUE)
        for v in self.cal.query('rid==9').drop_duplicates('id_pt').itertuples():
            col=BLUE if v.mis=='5c' else GREEN;ax.scatter(v.s_b_m/1000,v.z_b_m,c=col,marker='v');ax.scatter(v.s_pt_m/1000,v.H_dem_m,c=col)
        self.star(ax,self.k.loc[9].s_m/1000,self.k.loc[9].z_m,9);ax.set(xlabel='Distancia desde outlet (km)',ylabel='Elevación DEM (m)');self.save(fig,14,'Solo panel de perfil. Los dos paneles cartográficos están fuera del alcance.')

    def figura_15(self):
        q=self.means.merge(self.read('latitudes_calibracion.csv'),on=['rid','mis']);fig,axes=self.pair()
        for ax,m in zip(axes.flat,['G','P']):
            for mis,col,mark in [('5c',BLUE,'s'),('5e',GREEN,'o')]:
                c=q[q.mis.eq(mis)];ax.scatter(c.latitud_sur,c[f'K_{m}_m02_anio']*1e5,c=col,marker=mark,label='MIS '+mis)
                for v in c.itertuples():ax.annotate(f'R{v.rid:02d}',(v.latitud_sur,getattr(v,f'K_{m}_m02_anio')*1e5),xytext=(4,6),textcoords='offset points',fontsize=7)
            ax.set(title=f'K_{m}',xlabel='Latitud (°S)',ylabel='K (×10⁻⁵ m⁰·²/año)');ax.legend()
        self.save(fig,15)

    def figura_16(self):
        path=self.private/'comparacion_K_externa.csv'
        if not path.is_file():raise FileNotFoundError('comparacion_K_externa.csv')
        external=pd.read_csv(path);fig,ax=plt.subplots(figsize=(9,7),layout='constrained');ys=[];labels=[];j=0
        for name,lo,hi in [('Granitoides y metamórficas',1e-7,1e-6),('Volcaniclásticas',1e-5,1e-4),('Mudstones',1e-4,1e-2)]:ax.plot([lo,hi],[j,j],c=GRAY,lw=6);ys.append(j);labels.append(name+' · Stock y Montgomery (1999)');j+=1
        for v in external.itertuples():ax.errorbar(v.K,j,xerr=v.error_K,fmt='D',c=GREEN,capsize=3);ys.append(j);labels.append('Mataquito · '+v.sector+' · Long, com. pers.');j+=1
        for v in self.means.itertuples():
            ax.scatter(v.K_G_m02_anio,j-.1,c=BLUE,marker='o');ax.scatter(v.K_P_m02_anio,j+.1,c=GREEN,marker='s');ys.append(j);labels.append(f'R{v.rid:02d} · MIS {v.mis}');j+=1
        ax.set(xscale='log',yticks=ys,yticklabels=labels,xlabel='K (m⁰·²/año)');ax.invert_yaxis();self.save(fig,16,'Comparación externa suministrada localmente; valores no incluidos en el repositorio.')

    def figura_18(self):
        fig,axes=self.pair()
        for ax,mis in zip(axes.flat,['5c','5e']):
            x,v=self.timecurve(6,mis);self.curve(ax,x,self.p[6].elevacion_m,c=BLUE);self.star(ax,v.tau_G_ka,v.z_kp_m,6)
            tc=float(a.tiempo(self.c2,v.chi_origen_m,v.K_G_inferior_m02_anio,v.K_G_superior_m02_anio,self.c2));ax.scatter(tc,np.interp(self.c2,self.p[6].chi_m,self.p[6].elevacion_m),marker='D',c=GREEN,label='C2');ax.set(title='Donante MIS '+mis,xlabel='τ desde c (ka)',ylabel='Elevación DEM (m)');ax.legend()
        self.save(fig,18)

    def figura_19(self):
        q=a.calibrar_reloca();q=q[q.admisible];fig,axes=self.pair()
        for mis,col,dy in [('5c',BLUE,-.07),('5e',GREEN,.07)]:
            v=q[q.mis_donante.eq(mis)];y=v.frontera.str[1:].astype(int)+dy;axes[0,0].scatter(v.tau_kp_ka/1000,y,c=col,label='MIS '+mis);axes[0,1].scatter(v.RMSE_geometrico_ka,y,c=col)
        for ax in axes.flat:ax.set(yticks=[2,3,4],yticklabels=['C2','C3','C4'])
        axes[0,0].set_xlabel('Tiempo hasta KP (Ma)');axes[0,0].legend();axes[0,1].set_xlabel('RMSE geométrico B–Pₜ (ka)');self.save(fig,19)

    def figura_20(self):
        fig,axes=self.pair(2)
        for row,mis in enumerate(['5c','5e']):
            for ax,m in zip(axes[row],['G','P']):
                for v in self.t[self.t.mis.eq(mis)].itertuples():self.point(ax,v.latitud_sur,getattr(v,'tau_'+m+'_ka'),v.rid)
                ax.set(title=f'MIS {mis} · K_{m}',xlabel='Latitud (°S)',ylabel='τ desde c (ka)',ylim=(0,1950))
        self.save(fig,20)

    def figura_21(self):
        fig,axes=self.pair()
        for r,k in self.k.iterrows():self.point(axes[0,0],self.lat.loc[r],k.z_m,r)
        for v in self.t.query("mis=='5e'").itertuples():self.point(axes[0,1],v.tau_G_ka,v.z_kp_m,v.rid)
        axes[0,0].set(xlabel='Latitud (°S)',ylabel='Elevación KP (m)');axes[0,1].set(xlabel='τ con K_G–5e (ka)',ylabel='Elevación KP (m)');self.save(fig,21)

    def figura_22(self):
        data=self.read('sensibilidad_tiempos_altura.csv');fig,axes=self.pair()
        for ax,mis in zip(axes.flat,['5c','5e']):
            for y,r in enumerate(range(10,19)):
                q=data.query('rid==@r and mis==@mis');values=[q[q.fuente_altura.eq(source)].tau_G_ka.iloc[0]/1000 for source in ['DEM','JM']];ax.plot(values,[y,y],c=GRAY)
                for val,col,mark in zip(values,[BLUE,GREEN],['o','D']):ax.scatter(val,y,edgecolors=col,facecolors='white' if r in [11,14,17] else col,marker=mark)
            ax.set(yticks=range(9),yticklabels=[f'R{r:02d}' for r in range(10,19)],xlabel='τ al KP (Ma)',title='MIS '+mis);ax.invert_yaxis()
        self.save(fig,22)

    def figura_23(self):
        fig,axes=self.pair(3)
        for row,r in enumerate([2,3,4]):
            extra=self.read(f'R{r:02d}_completo.csv')
            for ax,col in zip(axes[row],['distancia_m','chi_m']):
                for p,ls,label in [(extra,'--','Completo'),(self.p[r],'-','Adoptado')]:ax.plot(p[col]/1000,p.elevacion_m,ls=ls,c=GRAY if ls=='--' else BLUE,label=label)
                ax.set(title=f'R{r:02d}',xlabel='Distancia (km)' if col=='distancia_m' else 'χ (km)',ylabel='Elevación DEM (m)');ax.legend(fontsize=7)
        self.save(fig,23)

    def figura_24(self):
        data=self.read('norte_donantes.csv');data=data[data.metodo.str.contains('PT')];fig,axes=plt.subplots(1,3,figsize=(10,4),sharey=True,layout='constrained')
        donors=list(data.donante.unique())
        for ax,r in zip(axes,[2,3,4]):
            q=data[data.rid.eq(r)]
            for j,donor in enumerate(donors):
                v=q[q.donante.eq(donor)].iloc[0];ax.scatter(j,v.tau_previo_ka,marker='*',c=RED,s=70);ax.scatter(j,v.tau_nuevo_ka,facecolors='white',edgecolors=BLUE)
            ax.axhline(125,c=GRAY,ls='--');ax.set(title=f'R{r:02d}',xticks=range(len(donors)),xticklabels=[s.replace(' MIS ','\n') for s in donors],xlabel='K del donante')
        axes[0].set_ylabel('τ del ensayo hasta KP (ka)');self.save(fig,24,'Escenarios del ensayo norte conservados; no sustituyen las selecciones finales.')

    def figura_26(self):
        fig,ax=plt.subplots(figsize=(9,4),layout='constrained');p=self.p[8];ax.plot(p.distancia_m/1000,p.elevacion_m,c=BLUE)
        cross=self.read('cruces_perfil_fallas.csv').query('rid==8 and falla=="Pelluhue"')
        for v in cross.itertuples():ax.axvline(v.s_cruce_m/1000,c=GRAY,ls=':')
        b=self.cal.query('rid==8').iloc[0];ax.scatter(b.s_b_m/1000,b.z_b_m,marker='v',c=BLUE,label='B');ax.scatter(p.distancia_m.iloc[0]/1000,p.elevacion_m.iloc[0],marker='s',c='black',label='Outlet SHP');self.star(ax,self.k.loc[8].s_m/1000,self.k.loc[8].z_m,8,label='KP');ax.set(xlabel='Distancia desde outlet SHP (km)',ylabel='Elevación DEM (m)');ax.legend();self.save(fig,26)

    def figura_28(self):
        fig,axes=self.pair()
        for ax,mis in zip(axes.flat,['5c','5e']):
            for v in self.rates.query('mis==@mis and ventana_solicitada_m==2000').itertuples():
                partial=v.fraccion_cubierta<.95;ax.scatter(self.lat.loc[v.rid],v.E_mm_anio,marker='D' if partial else 'o',facecolors='white' if partial or v.rid in [11,14,17] else BLUE,edgecolors=BLUE);ax.annotate(f'R{v.rid:02d}',(self.lat.loc[v.rid],v.E_mm_anio),xytext=(3,5),textcoords='offset points',fontsize=7)
            ax.set(title='MIS '+mis,xlabel='Latitud (°S)',ylabel='Incisión K_G ksn (mm/año)',ylim=(0,4.6))
        self.save(fig,28)

    def figura_29(self):
        path=self.private/'terrazas_comparacion.csv'
        if not path.is_file():raise FileNotFoundError('terrazas_comparacion.csv')
        jm=pd.read_csv(path);fig,axes=plt.subplots(2,1,figsize=(8,8),sharex=True,gridspec_kw={'height_ratios':[1,2]},layout='constrained');u=self.u.query("mis=='5e' and metodo=='G'")
        for ax,max_y in zip(axes,[3.3,1.1]):
            ax.errorbar(jm.latitud_sur,jm.U_mm_anio,yerr=jm.error_U_mm_anio,fmt='s',c=GREEN,alpha=.65,ms=4,capsize=2,label='Jara-Muñoz et al. (2015), 5e')
            for v in u.itertuples():
                ax.scatter(v.latitud_sur,v.U_KP_mm_anio,facecolors='white' if v.rid in [11,14,17] else BLUE,edgecolors=BLUE,zorder=5)
                if v.U_KP_mm_anio<max_y:ax.annotate(f'R{v.rid:02d}',(v.latitud_sur,v.U_KP_mm_anio),xytext=(3,5),textcoords='offset points',fontsize=7)
            ax.set(ylim=(min(0,(jm.U_mm_anio-jm.error_U_mm_anio).min())-.04,max_y),ylabel='Tasa o estimación (mm/año)')
        axes[0].legend(fontsize=8);axes[1].set_xlabel('Latitud (°S) · norte a sur');self.save(fig,29,'La serie costera se suministra localmente; U_KP usa desnivel y tiempo entre los mismos extremos.')

    def figura_30(self):
        fig,axes=plt.subplots(2,1,figsize=(8,6),layout='constrained');p=self.p[6]
        for ax,mask in zip(axes,[np.ones(len(p),dtype=bool),p.chi_m.between(self.cuts.chi_m.min()-500,self.k.loc[6].chi_m+500)]):
            ax.plot(p.chi_m[mask]/1000,p.elevacion_m[mask],c=BLUE)
            for v in self.cuts.itertuples():ax.axvline(v.chi_m/1000,c=BLUE if v.frontera=='C2' else GRAY,ls='--');ax.annotate(v.frontera,(v.chi_m/1000,v.z_m),xytext=(4,7),textcoords='offset points')
            self.star(ax,self.k.loc[6].chi_m/1000,self.k.loc[6].z_m,6);ax.set(xlabel='χ desde outlet (km)',ylabel='Elevación DEM (m)')
        self.save(fig,30)

    def figura_31(self):
        s=a.resumen_sur();fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
        for j,v in enumerate(s.itertuples()):ax.plot([v.minimo_ka/1000,v.maximo_ka/1000],[j,j],c=GRAY,lw=3);ax.scatter(v.promedio_ka/1000,j,c=BLUE);ax.scatter(v.mediana_ka/1000,j,marker='|',c=RED,s=130)
        ax.set(yticks=range(4),yticklabels=[f'K_{v.metodo} · MIS {v.mis}' for v in s.itertuples()],xlabel='τ de seis receptores reconocidos (Ma)');self.save(fig,31)

    def figura_33(self):
        manual=self.read('kp_referencias_manuales.csv').set_index('rid');base=self.t.query("mis=='5e' and rid in @a.SOUTH");alt=self.read('sensibilidad_tiempos_altura.csv');values=[[],[],[],[],[]]
        for v in base.itertuples():
            previous=float(a.tiempo(manual.loc[v.rid].chi_m_antes_revision,v.chi_origen_m,v.K_G_inferior_m02_anio));extended=float(a.tiempo(v.chi_kp_m,self.p[v.rid].chi_m.iloc[0],v.K_G_inferior_m02_anio))
            jm=alt.query('rid==@v.rid and mis=="5e" and fuente_altura=="JM"').tau_G_ka.iloc[0];c=self.t.query('rid==@v.rid and mis=="5c"').tau_G_ka.iloc[0]
            for bucket,z in zip(values,[100*(v.tau_G_ka/previous-1),100*(v.tau_G_ka/extended-1),100*(jm/v.tau_G_ka-1),100*(c/v.tau_G_ka-1),100*(v.tau_P_ka/v.tau_G_ka-1)]):bucket.append(z)
        fig,ax=plt.subplots(figsize=(9,4),layout='constrained')
        for j,z in enumerate(values):ax.plot([min(z),max(z)],[j,j],c=GRAY,lw=2);ax.scatter(z,[j]*len(z),c=BLUE)
        ax.axvline(0,c=GRAY,ls='--');ax.set(yticks=range(5),yticklabels=['KP manual → adoptado','Recorrido extendido → tramo','H DEM → publicada','MIS 5e → 5c','K_G → K_P'],xlabel='Cambio del tiempo (%)');ax.invert_yaxis();self.save(fig,33)

    def figura_34(self):
        fig,axes=plt.subplots(4,1,figsize=(11,11),layout='constrained',sharex=True,sharey=True)
        for ax,(lo,hi,mis) in zip(axes,[(1,9,'5c'),(1,9,'5e'),(10,18,'5c'),(10,18,'5e')]):
            for r in range(lo,hi+1):
                if self.t.query('rid==@r and mis==@mis').empty:continue
                x,v=self.timecurve(r,mis);self.curve(ax,x,self.p[r].elevacion_m,c=PALETTE[r-1],label=f'R{r:02d}');self.star(ax,v.tau_G_ka,v.z_kp_m,r)
            ax.set(title=f'R{lo:02d}–R{hi:02d} · MIS {mis}',xlabel='τ desde c (ka)',ylabel='Elevación DEM (m)');ax.legend(ncol=9,fontsize=7)
        self.save(fig,34)

    def figura_35(self):
        d=self.read('concavidad_exploratoria.csv');fig,axes=plt.subplots(2,1,figsize=(9,6),layout='constrained')
        for rid,q in d.groupby('ID'):axes[0].plot(q.theta,q.R2,label=f'R{int(rid):02d}',lw=1)
        best=d.loc[d.groupby('ID').R2.idxmax()];axes[0].axvline(.4,c='black',ls='--');axes[0].set(xlabel='θ ensayado',ylabel='R² z–χ');axes[0].legend(ncol=6,fontsize=7)
        axes[1].scatter(best.ID,best.theta,c=BLUE);axes[1].axhline(.4,c=GRAY,ls='--');axes[1].set(xticks=range(1,19),xlabel='Río',ylabel='θ con máximo R²');self.save(fig,35,'Diagnóstico guardado, previo a algunas geometrías adoptadas; no selecciona automáticamente θ.')

    def perfil_doble(self,n,r):
        fig,axes=self.pair();p=self.p[r]
        for ax,col in zip(axes.flat,['distancia_m','chi_m']):
            ax.plot((p[col]-p[col].iloc[0])/1000,p.elevacion_m,c=BLUE)
            if r in self.k.index:
                position=self.k.loc[r].s_m if col=='distancia_m' else self.k.loc[r].chi_m
                self.star(ax,(position-p[col].iloc[0])/1000,self.k.loc[r].z_m,r)
            ax.set(title=f'R{r:02d}'+(' · sin KP' if r==7 else ''),xlabel='Distancia desde outlet (km)' if col=='distancia_m' else 'χ desde outlet (km)',ylabel='Elevación DEM (m)')
        self.save(fig,n)

    def figura_39(self):
        d=self.read('JM2015_OD_alicuotas.csv');fig,ax=plt.subplots(figsize=(7,5),layout='constrained')
        for v in d.itertuples():
            selected=bool(v.seleccion_en_area);low=v.OD_porcentaje<=15;ax.scatter(v.OD_porcentaje,v.n_aliquotas,marker='s' if selected else 'o' if low else 'x',c=BLUE if selected else GREEN if low else GRAY);ax.annotate(v.muestra,(v.OD_porcentaje,v.n_aliquotas),xytext=(4,4),textcoords='offset points',fontsize=8)
        ax.axvline(15,c=GRAY,ls='--');ax.set(xlabel='Sobredispersión OD (%)',ylabel='Alícuotas aceptadas');self.save(fig,39,'Doce muestras de la tabla 1 publicada; sin datos del libro privado Ages_Final.')

    def ventanas(self,n,r):
        fig,axes=self.pair(2);p=self.p[r];k=self.k.loc[r];window=p.distancia_m.between(k.s_m-4100,k.s_m+350);axes[0,0].plot(p.distancia_m[window]/1000,p.elevacion_m[window],c=BLUE);self.star(axes[0,0],k.s_m/1000,k.z_m,r)
        rates=self.rates[self.rates.rid.eq(r)];q=rates[rates.mis.eq('5e')].sort_values('ventana_solicitada_m')
        for j,v in enumerate(q.itertuples()):
            mask=p.distancia_m.between(v.s_min_m-1e-7,v.s_max_m+1e-7);sub=p[mask];slope=v.ksn_m08/a.FACTOR;intercept=sub.elevacion_m.mean()-slope*sub.chi_m.mean();axes[0,1].plot(sub.chi_m/1000,slope*sub.chi_m+intercept,c=plt.colormaps['viridis'](j/7),label=f'{v.ventana_solicitada_m/1000:g} km')
        for mis,col in [('5c',BLUE),('5e',GREEN)]:
            q=rates[rates.mis.eq(mis)].sort_values('ventana_solicitada_m')
            if q.empty:continue
            axes[1,0].plot(q.ventana_solicitada_m/1000,q.E_mm_anio,c=col,label='MIS '+mis)
            for v in q.itertuples():axes[1,0].scatter(v.ventana_solicitada_m/1000,v.E_mm_anio,facecolors=col if v.fraccion_cubierta>=.95 else 'white',edgecolors=col)
            axes[1,1].plot(q.ventana_solicitada_m/1000,q.fraccion_cubierta*100,c=col,marker='o')
        axes[0,0].set(xlabel='Distancia del perfil (km)',ylabel='Elevación DEM (m)');axes[0,1].set(xlabel='χ (km)',ylabel='Elevación ajustada (m)');axes[0,1].legend(ncol=2,fontsize=7)
        axes[1,0].set(xlabel='Ventana solicitada (km)',ylabel='Incisión modelada (mm/año)');axes[1,0].legend();axes[1,1].set(xlabel='Ventana solicitada (km)',ylabel='Cobertura (%)',ylim=(0,105));axes[1,1].axhline(95,c=GRAY,ls='--');self.save(fig,n,'Pendientes guardadas; intercepto reconstruido sobre los mismos nodos del ajuste.')

    def figura_44(self):
        fig,axes=self.pair();im=None
        for ax,mis in zip(axes.flat,['5c','5e']):
            d=self.rates[self.rates.mis.eq(mis)];table=d.pivot(index='rid',columns='ventana_solicitada_m',values='cambio_vs_2km_pct')
            im=ax.pcolormesh(np.arange(len(table.columns)+1)-.5,np.arange(len(table)+1)-.5,table.to_numpy(),cmap='Blues',vmin=-100,vmax=100,shading='flat')
            ax.set(yticks=range(len(table)),yticklabels=[f'R{r:02d}'+(' ?' if r in [11,14,17] else '') for r in table.index],xticks=range(len(table.columns)),xticklabels=[f'{c/1000:g}' for c in table.columns],xlabel='Ventana solicitada (km)',title='MIS '+mis,ylim=(len(table)-.5,-.5));ax.grid(False)
        colorbar=fig.colorbar(im,ax=list(axes.flat),label='Cambio respecto de 2 km (%)',extend='both')
        colorbar.solids.set_rasterized(False)
        self.save(fig,44,'Escala saturada en ±100 %, con los valores completos conservados en CSV.')

    def figura_45(self):
        fig,axes=self.pair();p=self.p[6];v=self.t.query('rid==6 and mis=="5e"').iloc[0];axes[0,0].plot(p.chi_m/1000,p.elevacion_m,c=BLUE);axes[0,0].axvspan(0,v.chi_origen_m/1000,color='#e7edf2');self.star(axes[0,0],v.chi_kp_m/1000,v.z_kp_m,6)
        uniform=float(self.means.query('rid==6').K_G_m02_anio.iloc[0])
        for ki,ks,col,label in [(uniform,None,BLUE,'Uniforme'),(v.K_G_inferior_m02_anio,v.K_G_superior_m02_anio,GREEN,'C2')]:
            x=a.tiempo(p.chi_m,v.chi_origen_m,ki,ks,self.c2);self.curve(axes[0,1],x,p.elevacion_m,c=col,label=label);self.star(axes[0,1],float(a.tiempo(v.chi_kp_m,v.chi_origen_m,ki,ks,self.c2)),v.z_kp_m,6)
        axes[0,0].set(xlabel='χ desde outlet (km)',ylabel='Elevación DEM (m)');axes[0,1].set(xlabel='τ desde c (ka)',ylabel='Elevación DEM (m)');axes[0,1].legend();self.save(fig,45)

    def temporal_doble(self,n,r):
        fig,axes=self.pair()
        for ax,m in zip(axes.flat,['G','P']):
            for mis,col in [('5c',BLUE),('5e',GREEN)]:
                x,v=self.timecurve(r,mis,m);self.curve(ax,x,self.p[r].elevacion_m,c=col,label='MIS '+mis);self.star(ax,v[f'tau_{m}_ka'],v.z_kp_m,r)
            ax.set(title=f'R{r:02d} · K_{m}',xlabel='τ desde c (ka)',ylabel='Elevación DEM (m)');ax.legend()
        self.save(fig,n)

    def regresion(self,n,r):
        fits=self.read('regresiones_completo_3km.csv').set_index('rid').loc[r];p=self.p[r];k=self.k.loc[r];s=p.distancia_m.to_numpy();x=p.chi_m.to_numpy()/1000;fig,axes=self.pair();cut=float(fits.limite_inferior_ajuste_s_m);cutchi=np.interp(cut,p.distancia_m,p.chi_m)/1000
        for ax,prefix,title in zip(axes.flat,['completo','3km'],['Todo el tramo','3 km por lado']):
            ax.plot(x,p.elevacion_m,c=GRAY,lw=.8);ax.axvspan(x[0],cutchi,color='#e7edf2');ax.scatter(fits.chi_referencia_manual_m/1000,fits.z_referencia_manual_m,marker='+',c='black');self.star(ax,k.chi_m/1000,k.z_m,r)
            for side,col in [('inferior',BLUE),('superior',GREEN)]:
                slope=fits[f'{prefix}_{side}_b'];intercept=fits[f'{prefix}_{side}_a'];mask=(s>=fits[f'{prefix}_{side}_s_min_m']-1e-6)&(s<=fits[f'{prefix}_{side}_s_max_m']+1e-6)
                assert not np.any(s[mask]<cut-1e-5)
                ax.plot(x[mask],slope*p.chi_m[mask]+intercept,c=col,lw=1.7)
            inside=str(fits[f'{prefix}_cruce_inside']).lower()=='true';local=str(fits[f'{prefix}_cruce_local']).lower()=='true';chi=fits[f'{prefix}_cruce_chi']
            if inside and (local or prefix=='completo') and np.isfinite(chi) and chi>=cutchi*1000:ax.scatter(chi/1000,fits[f'{prefix}_cruce_z_profile'],marker='D',facecolors='white',edgecolors=GREEN)
            else:ax.text(.04,.06,'Cruce fuera del tramo evaluado',transform=ax.transAxes,fontsize=7)
            ax.set(title=f'R{r:02d} · {title}',xlabel='χ desde outlet (km)',ylabel='Elevación DEM (m)',xlim=(x.min()-.02*x.max(),x.max()*1.03),ylim=(p.elevacion_m.min()-15,p.elevacion_m.max()+15))
        self.save(fig,n,'Coeficientes OLS guardados de la revisión; cruce sobre el perfil, no cota del cruce de rectas.')

    def generar(self,seleccion=None):
        inventory=json.loads((DATA/'inventario_figuras.json').read_text(encoding='utf-8'))
        routes={36:(self.perfil_doble,1),37:(self.perfil_doble,7),38:(self.perfil_doble,17),41:(self.ventanas,3),42:(self.ventanas,6),43:(self.ventanas,8),46:(self.temporal_doble,14),47:(self.temporal_doble,9),48:(self.regresion,3),49:(self.regresion,8),50:(self.regresion,6),51:(self.regresion,9),52:(self.regresion,10),53:(self.regresion,14),54:(self.regresion,16),55:(self.regresion,18),56:(self.regresion,12)}
        for entry in inventory:
            n=entry['figura']
            if seleccion and n not in seleccion:continue
            if n in MAPS:self.record.append(dict(figura=n,estado='mapa_excluido'));continue
            try:
                if n==2:
                    p=self.private/'Jara_Munoz_2015_Fig12.png'
                    if not p.exists():raise FileNotFoundError(p.name)
                    target=self.out/'Figura_02_referencia_bibliografica.png';shutil.copy2(p,target);self.record.append(dict(figura=n,estado='reproduccion_bibliografica',archivos=[target.name],fuente='Jara-Muñoz et al. 2015, figura 12, p. 186; no calculada por este código.'))
                elif n in [4,5,6,7]:self.conceptual(n)
                elif n in routes:fn,r=routes[n];fn(n,r)
                else:getattr(self,f'figura_{n}')()
                print('OK',n,flush=True)
            except FileNotFoundError as e:
                if n not in [2,16,29]:raise
                self.record.append(dict(figura=n,estado='requiere_entrada_local',entrada=Path(str(e)).name));print('ENTRADA LOCAL',n,e,flush=True)
        (self.out/'cobertura.json').write_text(json.dumps(self.record,ensure_ascii=False,indent=2),encoding='utf-8')
        sources={p.relative_to(a.BASE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in a.DATA.rglob('*') if p.is_file()}
        for name in ['comparacion_K_externa.csv','terrazas_comparacion.csv','Jara_Munoz_2015_Fig12.png']:
            p=self.private/name
            if p.is_file():sources['entrada_local/'+name]=hashlib.sha256(p.read_bytes()).hexdigest()
        (self.out/'fuentes_graficos.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
        return self.record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--entradas-privadas',type=Path);p.add_argument('--figuras',type=int,nargs='+');p.add_argument('--formatos',nargs='+',choices=['png','svg','pdf'],default=['png','svg']);args=p.parse_args()
    output=args.output or a.BASE/'resultados'/('graficos_informe_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ'))
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42}):
        record=Graficos(output,args.entradas_privadas,args.formatos).generar(args.figuras)
    print(json.dumps({'carpeta':str(output),'generadas':sum(r['estado']=='generada' for r in record),'mapas_excluidos':sum(r['estado']=='mapa_excluido' for r in record),'entradas_pendientes':[r['figura'] for r in record if r['estado']=='requiere_entrada_local']},ensure_ascii=False))


if __name__=='__main__':main()
