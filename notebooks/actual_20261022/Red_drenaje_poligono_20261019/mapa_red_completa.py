"""Mapa de todas las ramas umbralizadas; visor Canvas autónomo, sin servidor adicional."""

import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
import base64,gzip,json,io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from PIL import Image
import figuras_red as fr

def crear_mapas(red,fallas,poly,dem,out,window=1000.,amin=.1):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    red=red[red.geom_type.eq('LineString')].copy()
    seg=np.array([np.array(g.coords)[[0,-1],:2] for g in red.geometry])
    good=np.isfinite(red.ksn)&red.ksn.gt(0)
    lo,hi=np.percentile(red.loc[good,'ksn'],[5,95]);lo=max(lo,.001);hi=max(hi,lo*1.01)
    norm=LogNorm(lo,hi,clip=True)
    bounds=poly.bounds;systems=fallas.nombre_base.unique().tolist()
    boxes={'Todo el polígono':list(bounds)}
    for name in systems:
        a,b,c,d=fallas[fallas.nombre_base.eq(name)].total_bounds
        boxes[name]=[max(bounds[0],a-2000),max(bounds[1],b-3000),min(bounds[2],c+2000),min(bounds[3],d+3000)]
    rgba=fr.hillshade(dem,bounds,maxside=3000)
    pngs=[]
    for name,box in boxes.items():
        fig,ax=plt.subplots(figsize=(12,11),layout='constrained')
        ax.imshow(rgba,extent=np.array(bounds)[[0,2,1,3]]/1000,origin='upper',alpha=.48)
        ax.add_collection(LineCollection(seg[~good]/1000,colors='#90969c',linewidths=.55,zorder=2))
        lc=LineCollection(seg[good]/1000,array=red.loc[good,'ksn'].to_numpy(),cmap='turbo',norm=norm,linewidths=1.05,zorder=3)
        ax.add_collection(lc)
        for system in systems:
            for geom in fallas[fallas.nombre_base.eq(system)].geometry:
                for line in (list(geom.geoms) if geom.geom_type=='MultiLineString' else [geom]):
                    xy=np.array(line.coords)/1000
                    ax.plot(xy[:,0],xy[:,1],color='white',lw=4.3,zorder=4)
                    ax.plot(xy[:,0],xy[:,1],color=fr.FALLA_COLORS.get(system,'black'),ls='--',lw=2,zorder=5)
        for geom in ([poly] if poly.geom_type=='Polygon' else poly.geoms):
            xy=np.array(geom.exterior.coords)/1000;ax.plot(xy[:,0],xy[:,1],color='#202020',lw=1.3,ls=':',zorder=6)
        ax.set_xlim(box[0]/1000,box[2]/1000);ax.set_ylim(box[1]/1000,box[3]/1000);ax.set_aspect('equal')
        ax.set_xlabel('Este UTM 18S (km)');ax.set_ylabel('Norte UTM 18S (km)');ax.grid(alpha=.16)
        ax.set_title(f'{name} · red de drenaje completa desde A ≥ {amin:g} km²\nksn en ventana de {window:g} m · θ = 0,40 · todas las ramas y órdenes',fontsize=14)
        fig.colorbar(lc,ax=ax,label='ksn (m⁰·⁸) · escala logarítmica · saturación visual P5–P95',shrink=.72,pad=.025)
        handles=[Line2D([],[],color=fr.FALLA_COLORS.get(s,'black'),ls='--',lw=2,label=s) for s in systems]
        handles.extend([Line2D([],[],color='#90969c',lw=2,label='Sin ksn válido'),Line2D([],[],color='#202020',ls=':',label='Polígono de estudio')])
        ax.legend(handles=handles,loc='upper left',fontsize=9,framealpha=.95)
        ax.annotate('N',xy=(.95,.91),xytext=(.95,.99),xycoords='axes fraction',ha='center',va='top',fontsize=13,
                    arrowprops={'arrowstyle':'<-','lw':1.7},zorder=10)
        width=(box[2]-box[0])/1000
        scale=10 if width>40 else 5 if width>15 else 2
        x=box[0]/1000+width*.04;y=box[1]/1000+(box[3]-box[1])/1000*.035
        ax.plot([x,x+scale],[y,y],color='white',lw=7,zorder=7);ax.plot([x,x+scale],[y,y],color='#111',lw=3,zorder=8)
        ax.text(x+scale/2,y+width*.012,f'{scale:g} km',ha='center',fontsize=10,bbox={'facecolor':'white','edgecolor':'none','alpha':.8},zorder=9)
        target=out/('Red_completa_ksn_'+name.replace(' ','_')+'.png')
        fig.savefig(target,dpi=180,facecolor='white');plt.close(fig);pngs.append(target)
    # Datos exactos de cada segmento, no un subconjunto de ríos ni un muestreo visual.
    # Coordenadas locales Float32 conservan precisión centimétrica en este polígono.
    xy=seg-np.array([bounds[0],bounds[1]])
    data=np.column_stack([xy.reshape(len(seg),4),red[['ksn','R2','area_km2','orden','span']].to_numpy()]).astype('<f4')
    payload=base64.b64encode(gzip.compress(data.tobytes(),compresslevel=6)).decode()
    buf=io.BytesIO();Image.fromarray(rgba).save(buf,format='PNG')
    bg='data:image/png;base64,'+base64.b64encode(buf.getvalue()).decode()
    colors=[plt.matplotlib.colors.to_hex(plt.get_cmap('turbo')(i/31)) for i in range(32)]
    faultlines=[]
    for name in systems:
        for geom in fallas[fallas.nombre_base.eq(name)].geometry:
            for line in (list(geom.geoms) if geom.geom_type=='MultiLineString' else [geom]):
                faultlines.append({'name':name,'color':fr.FALLA_COLORS.get(name,'black'),'xy':(np.array(line.coords)[:,:2]-np.array(bounds[:2])).tolist()})
    exterior=[(np.array(g.exterior.coords)-np.array(bounds[:2])).tolist() for g in ([poly] if poly.geom_type=='Polygon' else poly.geoms)]
    localboxes={k:[v[0]-bounds[0],v[1]-bounds[1],v[2]-bounds[0],v[3]-bounds[1]] for k,v in boxes.items()}
    cfg={'n':len(red),'stride':9,'width':bounds[2]-bounds[0],'height':bounds[3]-bounds[1],
         'x0':bounds[0],'y0':bounds[1],'lo':float(lo),'hi':float(hi),'window':window,'amin':amin,
         'colors':colors,'faults':faultlines,'poly':exterior,'boxes':localboxes}
    page=HTML.replace('__CFG__',json.dumps(cfg,ensure_ascii=False)).replace('__PAYLOAD__',payload).replace('__BG__',bg)
    htmlpath=out/'VER_RED_COMPLETA_KSN.html';htmlpath.write_text(page,encoding='utf8')
    return pngs,htmlpath,{'n_segmentos_dibujados':len(red),'n_ksn_validos':int(good.sum()),'color_min_P5':float(lo),'color_max_P95':float(hi),
                         'area_minima_km2':amin,'ventana_m':window,'sin_muestreo_geometria':True}

HTML=r'''<!doctype html><html lang="es"><meta charset="utf-8"><title>Red completa del polígono · ksn</title>
<style>*{box-sizing:border-box}body{margin:0;font:15px system-ui;color:#173348;background:#eef2f5}header{padding:15px 22px;background:white;border-bottom:1px solid #cbd5df}h1{font-size:22px;margin:0 0 7px}p{margin:6px 0}.toolbar{display:flex;gap:15px;align-items:center;flex-wrap:wrap;margin-top:12px}select,button{font:inherit;padding:7px;border:1px solid #abb9c8;border-radius:5px;background:white}button{cursor:pointer}.map{position:relative;margin:12px;background:white;border:1px solid #cbd5df;height:calc(100vh - 205px);min-height:500px}canvas{display:block;width:100%;height:100%;touch-action:none;cursor:grab}.legend{position:absolute;right:15px;top:15px;background:#fffffff2;border:1px solid #cad2dc;padding:12px;border-radius:5px;min-width:230px;pointer-events:none}.gradient{height:16px;border:1px solid #777;margin:8px 0}.ticks{display:flex;justify-content:space-between;font-size:12px}.fault{margin-top:8px;font-size:13px}.tip{display:none;position:absolute;pointer-events:none;background:#ffffffef;border:1px solid #8898a8;border-radius:5px;padding:9px;box-shadow:0 3px 10px #0002;z-index:10;font-size:13px;max-width:285px}.status{position:absolute;bottom:15px;left:15px;background:#fffffff2;padding:7px;border-radius:4px;pointer-events:none;font-size:13px}footer{padding:0 22px 15px;font-size:13px}</style>
<header><h1>Red de drenaje completa dentro de PoligonoFallas</h1><p id="intro">Cargando todos los segmentos…</p>
<div class="toolbar"><label>Vista <select id="view"></select></label><label>Escala de color <select id="mode"><option value="log">Logarítmica</option><option value="linear">Lineal</option></select></label>
<label>Área mínima visible <select id="area"></select> km²</label><label><input id="relief" type="checkbox" checked> Hillshade</label><label><input id="faults" type="checkbox" checked> Fallas</label><label><input id="band" type="checkbox"> Resaltar ±1 km de las fallas</label><button id="reset">Ver todo</button><button id="save">Guardar vista PNG</button></div></header>
<div class="map" id="container"><canvas id="map"></canvas><div class="legend"><b>ksn (m⁰·⁸)</b><div id="gradient" class="gradient"></div><div class="ticks" id="ticks"></div><small>P5–P95; valores extremos conservados</small><div id="faultlegend"></div><div class="fault">Gris: sin ajuste válido<br>Línea punteada negra: polígono</div></div><div class="tip" id="tip"></div><div class="status" id="status">Preparando la red…</div></div>
<footer>Rueda: zoom · arrastrar: desplazar · pasar el cursor: valor de cada segmento. Todos los órdenes de cauce están incluidos desde el umbral indicado. La ventana se cambia en el notebook. El área aportante procede del DEM completo; las líneas se recortan al polígono.</footer>
<script>
const C=__CFG__, B64='__PAYLOAD__', bg=new Image();bg.src='__BG__';
const $=id=>document.getElementById(id),canvas=$('map'),ctx=canvas.getContext('2d'),box=$('container');
let D,paths=[],cam,screenW,screenH,scale=1,drag=null,moved=false,ready=false,frame=0,shown=0;
const buckets=new Map(),bin=400;
function fit(b){let w=b[2]-b[0],h=b[3]-b[1];cam={x:(b[0]+b[2])/2,y:(b[1]+b[3])/2,width:Math.max(w*1.04,h*screenW/screenH*1.04)};draw();}
function transform(){scale=screenW/cam.width;ctx.setTransform(scale*devicePixelRatio,0,0,-scale*devicePixelRatio,(screenW/2-cam.x*scale)*devicePixelRatio,(screenH/2+cam.y*scale)*devicePixelRatio);}
function world(px,py){return [cam.x+(px-screenW/2)/scale,cam.y-(py-screenH/2)/scale];}
function build(){paths=Array.from({length:33},()=>new Path2D());shown=0;let area=Number($('area').value),lg=$('mode').value==='log';
 for(let i=0;i<C.n;i++){let j=i*9;if(D[j+6]+1e-6<area)continue;shown++;let v=D[j+4],k=32;
 if(Number.isFinite(v)&&v>0){let t=lg?(Math.log(v)-Math.log(C.lo))/(Math.log(C.hi)-Math.log(C.lo)):(v-C.lo)/(C.hi-C.lo);k=Math.max(0,Math.min(31,Math.floor(t*32)));}
 paths[k].moveTo(D[j],D[j+1]);paths[k].lineTo(D[j+2],D[j+3]);}
 let ticks=Array.from({length:4},(_,i)=>lg?C.lo*Math.pow(C.hi/C.lo,i/3):C.lo+(C.hi-C.lo)*i/3);
 $('ticks').innerHTML=ticks.map(v=>'<span>'+v.toPrecision(3)+'</span>').join('');draw();}
function line(xy){ctx.beginPath();xy.forEach((p,i)=>i?ctx.lineTo(...p):ctx.moveTo(...p));}
function draw(){if(!ready||!cam)return;if(frame)return;frame=requestAnimationFrame(()=>{frame=0;ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);ctx.fillStyle='#f9fbfc';ctx.fillRect(0,0,screenW,screenH);transform();
 if($('relief').checked&&bg.complete){ctx.save();ctx.globalAlpha=.55;ctx.translate(0,C.height);ctx.scale(1,-1);ctx.drawImage(bg,0,0,C.width,C.height);ctx.restore();}
 if($('band').checked){ctx.lineWidth=2000;ctx.lineCap='round';ctx.strokeStyle='#6a22771a';C.faults.forEach(f=>{line(f.xy);ctx.stroke();});ctx.lineCap='butt';}
 ctx.lineWidth=1.65/scale;for(let i=0;i<33;i++){ctx.strokeStyle=i===32?'#999fa5':C.colors[i];ctx.stroke(paths[i]);}
 if($('faults').checked){C.faults.forEach(f=>{line(f.xy);ctx.lineWidth=6/scale;ctx.strokeStyle='white';ctx.setLineDash([]);ctx.stroke();ctx.lineWidth=2.5/scale;ctx.strokeStyle=f.color;ctx.setLineDash([9/scale,5/scale]);ctx.stroke();});}
 ctx.lineWidth=1.8/scale;ctx.strokeStyle='#20272e';ctx.setLineDash([3/scale,5/scale]);C.poly.forEach(p=>{line(p);ctx.stroke();});ctx.setLineDash([]);
 ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);ctx.fillStyle='#173348';ctx.font='bold 16px system-ui';ctx.fillText('N ↑',screenW-55,screenH-45);
 let raw=cam.width*.13,pow=10**Math.floor(Math.log10(raw)),s=[1,2,5,10].map(v=>v*pow).filter(v=>v<=raw).pop()||pow;
 ctx.strokeStyle='#111';ctx.lineWidth=4;ctx.beginPath();ctx.moveTo(28,screenH-65);ctx.lineTo(28+s*scale,screenH-65);ctx.stroke();ctx.font='13px system-ui';ctx.fillText((s/1000).toLocaleString('es-CL')+' km',28,screenH-75);
 $('status').textContent=shown.toLocaleString('es-CL')+' segmentos dibujados · ventana '+C.window+' m · EPSG:32718';});}
function resize(){screenW=box.clientWidth;screenH=box.clientHeight;canvas.width=screenW*devicePixelRatio;canvas.height=screenH*devicePixelRatio;if(!cam)fit(C.boxes['Todo el polígono']);else draw();}
function eventXY(e){let r=canvas.getBoundingClientRect();return [e.clientX-r.left,e.clientY-r.top];}
canvas.addEventListener('wheel',e=>{e.preventDefault();if(!ready)return;let [px,py]=eventXY(e),p=world(px,py),factor=Math.exp(Math.max(-.5,Math.min(.5,e.deltaY*.001)));cam.width=Math.max(60,Math.min(C.width*5,cam.width*factor));scale=screenW/cam.width;cam.x=p[0]-(px-screenW/2)/scale;cam.y=p[1]+(py-screenH/2)/scale;$('tip').style.display='none';draw();},{passive:false});
canvas.onpointerdown=e=>{if(!ready)return;drag={p:eventXY(e),x:cam.x,y:cam.y};moved=false;canvas.setPointerCapture(e.pointerId);canvas.style.cursor='grabbing';};
canvas.onpointerup=e=>{drag=null;canvas.style.cursor='grab';};
canvas.onpointerleave=()=>{$('tip').style.display='none';};
canvas.onpointermove=e=>{if(!ready)return;let [px,py]=eventXY(e);if(drag){cam.x=drag.x-(px-drag.p[0])/scale;cam.y=drag.y+(py-drag.p[1])/scale;moved=true;$('tip').style.display='none';draw();return;}
 let [x,y]=world(px,py),ix=Math.floor(x/bin),iy=Math.floor(y/bin),rad=9/scale,r=Math.max(1,Math.ceil(rad/bin)),best=rad*rad,hit=-1,area=+$('area').value;
 for(let a=ix-r;a<=ix+r;a++)for(let b=iy-r;b<=iy+r;b++){let ids=buckets.get(a+','+b)||[];for(let i of ids){let j=i*9;if(D[j+6]+1e-6<area)continue;let dx=D[j+2]-D[j],dy=D[j+3]-D[j+1],t=Math.max(0,Math.min(1,((x-D[j])*dx+(y-D[j+1])*dy)/(dx*dx+dy*dy||1))),dd=(x-D[j]-t*dx)**2+(y-D[j+1]-t*dy)**2;if(dd<best){best=dd;hit=i;}}}
 let tip=$('tip');if(hit<0){tip.style.display='none';return;}let j=hit*9,v=D[j+4];tip.innerHTML='<b>ksn '+(Number.isFinite(v)?v.toFixed(2)+' m⁰·⁸':'sin ajuste válido')+'</b><br>Área '+D[j+6].toFixed(3)+' km² · orden '+Math.round(D[j+7])+'<br>R² '+(Number.isFinite(D[j+5])?D[j+5].toFixed(3):'—')+'<br>Ventana efectiva '+Math.round(D[j+8])+' m<br>Segmento '+(hit+1)+'<br>UTM '+Math.round(C.x0+x)+', '+Math.round(C.y0+y);tip.style.left=Math.min(px+14,screenW-285)+'px';tip.style.top=Math.min(py+14,screenH-180)+'px';tip.style.display='block';};
$('view').innerHTML=Object.keys(C.boxes).map(n=>'<option>'+n+'</option>').join('');$('view').onchange=()=>fit(C.boxes[$('view').value]);
$('area').innerHTML=[C.amin,.25,.5,1,2,5].filter((v,i,a)=>v>=C.amin&&a.indexOf(v)===i).map(v=>'<option>'+v+'</option>').join('');$('area').onchange=build;$('mode').onchange=build;
['relief','faults','band'].forEach(id=>$(id).onchange=draw);$('reset').onclick=()=>{$('view').value='Todo el polígono';fit(C.boxes['Todo el polígono']);};
$('save').onclick=()=>{let a=document.createElement('a');a.download='Red_completa_ksn_vista.png';a.href=canvas.toDataURL('image/png');a.click();};
$('gradient').style.background='linear-gradient(to right,'+C.colors.join(',')+')';$('faultlegend').innerHTML=[...new Set(C.faults.map(f=>f.name))].map(n=>'<div class="fault" style="color:'+C.faults.find(f=>f.name===n).color+'">━ ━ '+n+'</div>').join('');
bg.onload=draw;new ResizeObserver(resize).observe(box);
(async()=>{try{let bytes=Uint8Array.from(atob(B64),c=>c.charCodeAt(0)),buf=await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();D=new Float32Array(buf);if(D.length!==C.n*9)throw Error('Tamaño de red inválido');
 for(let i=0;i<C.n;i++){let j=i*9,key=Math.floor((D[j]+D[j+2])*.5/bin)+','+Math.floor((D[j+1]+D[j+3])*.5/bin);if(!buckets.has(key))buckets.set(key,[]);buckets.get(key).push(i);}
 ready=true;$('intro').textContent='Todas las ramas con A ≥ '+C.amin+' km² · '+C.n.toLocaleString('es-CL')+' segmentos · ksn con ventana '+C.window+' m · θ=0,40';build();resize();}catch(e){$('status').textContent='No se pudo abrir el visor: '+e.message+' · utiliza las figuras PNG del notebook.';console.error(e);}})();
</script></html>'''
