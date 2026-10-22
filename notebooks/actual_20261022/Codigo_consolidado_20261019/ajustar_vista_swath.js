/* Recorte y estadísticas sobre el DEM nativo; EV sólo cambia proporciones. */
globalThis.SwathVista = (() => {
  'use strict';
  const cache=new WeakMap();
  const safe=v=>Number.isFinite(v)?v:null;
  function decode(s){const b=atob(s),u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);return new Float32Array(u.buffer);}
  function native(r){if(!cache.has(r))cache.set(r,{z:decode(r.z_b64),s:decode(r.pendiente_b64)});return cache.get(r);}
  function sample(r,a,x,y){
    const c=(x-r.x0)/r.dx,h=(y-r.y0)/r.dy,j=Math.floor(c),i=Math.floor(h);
    if(j<0||i<0||j+1>=r.nx||i+1>=r.ny)return NaN;
    const fx=c-j,fy=h-i,k=i*r.nx+j,v=[a[k],a[k+1],a[k+r.nx],a[k+r.nx+1]];
    if(!v.every(Number.isFinite))return NaN;
    return (v[0]*(1-fx)+v[1]*fx)*(1-fy)+(v[2]*(1-fx)+v[3]*fx)*fy;
  }
  function axis(size,step){let n=Math.max(3,Math.ceil(size/step)+1);if(n%2===0)n++;return Array.from({length:n},(_,i)=>-size/2+i*size/(n-1));}
  const percentile=(v,q)=>{if(!v.length)return null;const pos=(v.length-1)*q,i=Math.floor(pos),f=pos-i;return v[i]*(1-f)+v[Math.min(i+1,v.length-1)]*f;};
  function stats(v){const a=v.filter(Number.isFinite).sort((a,b)=>a-b);return {z_min_m:a.length?a[0]:null,z_max_m:a.length?a[a.length-1]:null,z_mediana_m:percentile(a,.5),z_p90_m:percentile(a,.9),coverage:a.length/v.length};}
  function adjust(view,params){
    const r=view.raster,c=view.cfg,{ancho:w,largo:l,ev}=params;
    if(![w,l,ev].every(Number.isFinite)||w<100||l<100||w>r.max_ancho_m||l>r.max_largo_m||ev<.1||ev>30)
      throw new Error('Ancho y largo: 100–'+r.max_ancho_m+' m. Exageración vertical: 0,1–30.');
    const a=native(r),east=axis(w,r.resolucion_m),north=axis(l,r.resolucion_m);
    const sampleZ=(x,y)=>sample(r,a.z,c.x_B+x,c.y_B+y);
    const full=north.map(y=>east.map(x=>sampleZ(x,y)));
    const ns=full.map(stats),ew=east.map((_,j)=>stats(full.map(row=>row[j])));
    const midx=(east.length-1)/2,midy=(north.length-1)/2;
    // Estadísticas a resolución nativa; superficie limitada para mantener la interacción.
    const factor=Math.max(1,Math.sqrt(east.length*north.length/90000));
    const ex=factor===1?east:axis(w,r.resolucion_m*factor),ny=factor===1?north:axis(l,r.resolucion_m*factor);
    const z=factor===1?full:ny.map(y=>ex.map(x=>sampleZ(x,y)));
    const slope=ny.map(y=>ex.map(x=>safe(sample(r,a.s,c.x_B+x,c.y_B+y))));
    const xx=ny.map(()=>ex.slice()),yy=ny.map(y=>ex.map(()=>y));
    let lo=Infinity,hi=-Infinity,valid=0;
    for(const row of full)for(const v of row)if(Number.isFinite(v)){lo=Math.min(lo,v);hi=Math.max(hi,v);valid++;}
    if(!valid)throw new Error('El recorte no contiene elevaciones válidas.');
    const data=structuredClone(view.fig.data),layout=structuredClone(view.fig.layout),zB=safe(sampleZ(0,0));
    const jmCounts={'5e':0,'5c':0};
    for(const t of data){
      const m=t.meta||{};
      if(m.capa==='DEM'){
        const colored=params.pendiente!==false;
        Object.assign(t,{x:xx,y:yy,z:z.map(row=>row.map(safe)),surfacecolor:colored?slope:z.map(row=>row.map(safe))});
        if(!colored){t.colorscale=[[0,'#656d77'],[.5,'#bec4ca'],[1,'#f1f3f5']];t.cmin=lo;t.cmax=hi;t.colorbar={...(t.colorbar||{}),title:{text:'Elevación DEM (m)'}};}
      }
      if(m.capa==='B_inicial'){t.x=[0];t.y=[0];t.z=[zB];t.text=['B inicial'];}
      if(m.capa==='B_perfil'){t.x=[0];t.y=[zB];}
      if(m.capa==='rio'&&c.river){
        const p=c.river,ok=p.x.map((v,i)=>Number.isFinite(v)&&Number.isFinite(p.y[i])&&Math.abs(v)<=w/2&&Math.abs(p.y[i])<=l/2);
        for(const k of ['x','y','z'])t[k]=p[k].map((v,i)=>ok[i]?v:null);
      }
      if(m.capa==='terrazas_JM'){
        const pts=(c.JM_catalogo||[]).filter(p=>p.mis_vis===m.mis&&Math.abs(p.este_utm18s_m-c.x_B)<=w/2&&Math.abs(p.norte_utm18s_m-c.y_B)<=l/2)
          .map(p=>({...p,x:p.este_utm18s_m-c.x_B,y:p.norte_utm18s_m-c.y_B,z:sampleZ(p.este_utm18s_m-c.x_B,p.norte_utm18s_m-c.y_B)})).filter(p=>Number.isFinite(p.z));
        t.x=pts.map(p=>p.x);t.y=pts.map(p=>p.y);t.z=pts.map(p=>p.z+3);t.text=pts.map(p=>p.terrace_id);
        t.hovertext=pts.map(p=>`${p.terrace_id} · MIS ${p.mis_vis}<br>Altura publicada: ${p.elevacion_m.toFixed(2)} ± ${p.error_elevacion_m.toFixed(2)} m<br>Edad MIS: ${p.edad_referencia_ka} ± ${p.incertidumbre_edad_ka} ka<br>Este UTM: ${p.este_utm18s_m.toFixed(1)} m<br>Norte UTM: ${p.norte_utm18s_m.toFixed(1)} m<br>Terreno DEM en XY: ${p.z.toFixed(2)} m<br>Marcador sobre DEM (+3 m visuales)`);
        t.meta.n_puntos=pts.length;jmCounts[m.mis]=pts.length;
      }
      if(m.capa==='eje_NS'){t.x=north.map(()=>0);t.y=north;t.z=full.map(row=>safe(row[midx]));}
      if(m.capa==='eje_EW'){t.x=east;t.y=east.map(()=>0);t.z=full[midy].map(safe);}
      if(m.capa==='perfil_swath'){const arr=m.direccion==='NS'?ns:ew;t.x=m.direccion==='NS'?north:east;t.y=arr.map(v=>v[m.variable]);t.customdata=arr.map(v=>v.coverage);}
      if(m.capa==='perfil_central'){t.x=m.direccion==='NS'?north:east;t.y=m.direccion==='NS'?full.map(row=>safe(row[midx])):full[midy].map(safe);}
      if(m.capa==='plano_H'){t.x=[[-w/2,w/2],[-w/2,w/2]];t.y=[[-l/2,-l/2],[l/2,l/2]];t.z=[[m.H_m,m.H_m],[m.H_m,m.H_m]];}
    }
    const span=Math.max(1,hi-lo),scale=Math.max(w,l,span*ev);
    Object.assign(layout.scene,{aspectmode:'manual',aspectratio:{x:w/scale,y:l/scale,z:span*ev/scale}});
    layout.scene.xaxis.range=[-w/2,w/2];layout.scene.xaxis.autorange=false;
    layout.scene.yaxis.range=[-l/2,l/2];layout.scene.yaxis.autorange=false;
    layout.title.text=`${c.rio} · ${c.caso_B}<br><sup>Bandas de relieve: ${w} m O–E × ${l} m S–N · ${params.pendiente!==false?'azul = menor pendiente; rojo = mayor':'relieve gris, sin color de pendiente'} · EV=${ev} · alturas DEM</sup>`;
    if(layout.annotations?.[1])layout.annotations[1].text=`Banda N–S · ancho O–E: ${w} m`;
    if(layout.annotations?.[2])layout.annotations[2].text=`Banda O–E · largo N–S: ${l} m`;
    const cfg={...c,east_m:ex,north_m:ny,z_m:z.map(row=>row.map(safe)),pendiente_grados:slope,ancho_ew_m:w,largo_ns_m:l,exageracion_vertical:ev,JM_en_3D:jmCounts};
    const contexto=view.contexto?structuredClone(view.contexto):null;
    if(contexto){
      Object.assign(contexto.layout.shapes[0],{x0:-w/2,x1:w/2,y0:-l/2,y1:l/2});
      for(const t of contexto.data)if(t.meta?.capa==='terrazas_JM')
        t.hovertext=t.hovertext.map((h,i)=>h.replace(/(?:Dentro|Fuera) del swath 3D/,Math.abs(t.x[i])<=w/2&&Math.abs(t.y[i])<=l/2?'Dentro del swath 3D':'Fuera del swath 3D'));
    }
    return {fig:{data,layout},cfg,contexto,info:`Cobertura DEM: ${(100*valid/(east.length*north.length)).toFixed(1)} % · reserva: ${r.max_ancho_m/1000} × ${r.max_largo_m/1000} km · muestreo ≈${r.resolucion_m} m`};
  }
  return {ajustarVista:adjust};
})();
