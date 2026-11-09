(async function () {
  'use strict';
  const compressed=Uint8Array.from(atob(document.getElementById('datos-visor').textContent),c=>c.charCodeAt(0));
  const dataset = JSON.parse(await new Response(new Blob([compressed]).stream().pipeThrough(new DecompressionStream('gzip'))).text());
  const cases = dataset.casos;
  const gd = document.getElementById('swath-map');
  const contextGD=document.getElementById('visor-contexto');
  const contextWrap=document.getElementById('visor-contexto-wrap');
  const terraInfo=document.getElementById('visor-terrazas-info');
  const width=document.getElementById('visor-ancho'),length=document.getElementById('visor-largo'),ev=document.getElementById('visor-ev');
  const apply=document.getElementById('visor-aplicar'),gridInfo=document.getElementById('visor-malla-info');
  const parameters=()=>({ancho:Number(width.value),largo:Number(length.value),ev:Number(ev.value),pendiente:document.getElementById('visor-color-pendiente').checked});
  const terraChecks={'5e':document.getElementById('visor-mis5e'),'5c':document.getElementById('visor-mis5c')};
  for(const mis of ['5e','5c'])terraChecks[mis].checked=dataset.mostrarTerrazas?.[mis]!==false;
  const rio = document.getElementById('visor-rio'), escenario = document.getElementById('visor-escenario');
  const message = document.getElementById('visor-status'), total = document.getElementById('visor-total');
  let active = null, busy = false, globalHandle = null, globalQueue = Promise.resolve();
  const byId = new Map(cases.map(c => [c.cfg.caso_B,c]));
  const remembered = new Map();
  const keyOf = cfg => 'swath_puntos_v1:' + cfg.caso_B + ':' + cfg.x_B.toFixed(3) + ':' + cfg.y_B.toFixed(3);
  const columns = ['id_punto','caso_B','rid','rio','mis','tipo','etiqueta','nota','x_utm18s_m','y_utm18s_m',
    'elevacion_dem_m','pendiente_grados','crs','fecha_utc','dem_fuente','metodo_posicion'];
  const quote = v => '"' + String(v == null ? '' : v).replace(/"/g,'""') + '"';
  function saved(c) {
    try {
      const v = JSON.parse(localStorage.getItem(keyOf(c.cfg)) || 'null');
      return v && Array.isArray(v.points) ? v.points.filter(p => p.caso_B===c.cfg.caso_B &&
        [p.x_utm18s_m,p.y_utm18s_m,p.elevacion_dem_m,p.pendiente_grados].every(Number.isFinite))
          .map(p => ({...p,rio:c.cfg.rio,rid:c.cfg.rid,mis:c.cfg.mis})) : [];
    } catch (_) { return []; }
  }
  function allPoints() { return cases.flatMap(saved); }
  function csv(rows) { return '\ufeff'+columns.join(',')+'\r\n'+rows.map(p => columns.map(c => quote(p[c])).join(',')).join('\r\n')+(rows.length?'\r\n':''); }
  function counts() {
    const n=allPoints().length;
    total.textContent=n+' punto(s) entre todos los ríos';
    return n;
  }
  function writeGlobal() {
    if (!globalHandle) return;
    const text=csv(allPoints()),handle=globalHandle;
    globalQueue=globalQueue.then(async()=>{
      const w=await handle.createWritable();
      try{await w.write(text);await w.close();}
      catch(e){try{await w.abort();}catch(_){}throw e;}
      message.textContent='CSV de todos los ríos actualizado: '+handle.name;
    }).catch(e=>{globalHandle=null;message.textContent='No se pudo escribir el CSV global; usa «Descargar todos». Los puntos se conservan.';});
  }
  function option(value,text){const o=document.createElement('option');o.value=String(value);o.textContent=text;return o;}
  const rivers=new Map(cases.map(c=>[String(c.cfg.rid),c.cfg.rio]));
  for(const [id,name] of [...rivers].sort((a,b)=>Number(a[0])-Number(b[0])))rio.appendChild(option(id,'R'+id.padStart(2,'0')+' · '+name));
  function fillCases(id,desired){
    escenario.replaceChildren();
    const group=cases.filter(c=>String(c.cfg.rid)===String(id));
    group.forEach(c=>escenario.appendChild(option(c.cfg.caso_B,c.etiqueta)));
    escenario.value=group.some(c=>c.cfg.caso_B===desired)?desired:group[0].cfg.caso_B;
  }
  function lock(on){busy=on;rio.disabled=on;escenario.disabled=on;apply.disabled=on;gd.style.pointerEvents=on?'none':'';}
  async function applyTerraces(){
    const jobs=[];
    for(const div of [gd,contextGD]){
      if(!div.data)continue;
      div.data.forEach((tr,i)=>{
        if(tr.meta?.capa==='terrazas_JM')jobs.push(Plotly.restyle(div,{visible:terraChecks[tr.meta.mis].checked?true:'legendonly'},[i]));
      });
    }
    await Promise.all(jobs);
  }
  function terraceInfo(view){
    const a=view.cfg.JM_en_3D||{'5e':0,'5c':0},b=view.cfg.JM_en_contexto;
    terraInfo.textContent='En 3D: '+a['5e']+' de 5e y '+a['5c']+' de 5c'+
      (b?' · en contexto: '+b['5e']+' de 5e y '+b['5c']+' de 5c':'');
  }
  for(const mis of ['5e','5c'])terraChecks[mis].addEventListener('change',applyTerraces);
  async function change(id){
    if(busy || !byId.has(id))return;
    lock(true);message.textContent='Cargando relieve y puntos de '+byId.get(id).cfg.rio+'…';
    try{
      const next=byId.get(id);
      const rendered=SwathVista.ajustarVista(next,parameters());
      if(active?.cfg.caso_B===id&&gd.layout?.scene?.camera)rendered.fig.layout.scene.camera=structuredClone(gd.layout.scene.camera);
      if(gd.swathPuntos)await gd.swathPuntos.prepararCambio();
      await Plotly.react(gd,rendered.fig.data,rendered.fig.layout,dataset.plotConfig);
      contextWrap.hidden=!rendered.contexto;
      if(rendered.contexto)await Plotly.react(contextGD,rendered.contexto.data,rendered.contexto.layout,dataset.plotConfig);
      await applyTerraces();terraceInfo(rendered);gridInfo.textContent=rendered.info;
      await gd.swathPuntos.cambiarCaso(rendered.cfg);
      active=next;remembered.set(String(next.cfg.rid),id);
      try{localStorage.setItem('swath_visor_caso_activo_v1',id);}catch(_){}
      message.textContent=next.cfg.rio+' · '+next.etiqueta+' · los puntos se guardan en este escenario.';
      counts();
    }catch(e){message.textContent='No se pudo aplicar la vista: '+e.message;}finally{lock(false);}
  }
  apply.addEventListener('click',()=>change(escenario.value));
  for(const input of [width,length,ev])input.addEventListener('keydown',e=>{if(e.key==='Enter')change(escenario.value);});
  rio.addEventListener('change',()=>{fillCases(rio.value,remembered.get(rio.value));change(escenario.value);});
  escenario.addEventListener('change',()=>change(escenario.value));
  document.getElementById('visor-descargar-todos').addEventListener('click',()=>{
    const rows=allPoints(),url=URL.createObjectURL(new Blob([csv(rows)],{type:'text/csv;charset=utf-8;'}));
    const a=document.createElement('a');a.href=url;a.download='Puntos_B_todos_los_rios.csv';document.body.appendChild(a);a.click();a.remove();
    setTimeout(()=>URL.revokeObjectURL(url),15000);message.textContent=rows.length+' punto(s) exportado(s) con río, MIS y escenario.';
  });
  document.getElementById('visor-vincular-global').addEventListener('click',async()=>{
    if(!window.showSaveFilePicker){message.textContent='Este navegador no permite escritura directa. Usa «Descargar todos». ';return;}
    try{
      globalHandle=await window.showSaveFilePicker({suggestedName:'Puntos_B_todos_los_rios.csv',
        types:[{description:'Puntos de todos los ríos',accept:{'text/csv':['.csv']}}]});
      writeGlobal();
    }catch(e){if(e.name!=='AbortError')message.textContent='No se pudo vincular el CSV global; usa «Descargar todos».';}
  });
  window.addEventListener('swath-puntos-cambiaron',()=>{counts();writeGlobal();});

  // CSV estándar: incluye comas y saltos de línea dentro de campos entrecomillados.
  function parseCSV(text){
    text=text.replace(/^\ufeff/,'');let rows=[],row=[],field='',quoted=false;
    for(let i=0;i<text.length;i++){
      const c=text[i];
      if(c==='"'){if(quoted&&text[i+1]==='"'){field+='"';i++;}else quoted=!quoted;}
      else if(c===','&&!quoted){row.push(field);field='';}
      else if((c==='\n'||c==='\r')&&!quoted){if(c==='\r'&&text[i+1]==='\n')i++;row.push(field);if(row.some(v=>v!==''))rows.push(row);row=[];field='';}
      else field+=c;
    }
    if(quoted)throw new Error('CSV con comillas sin cerrar');
    if(field!==''||row.length){row.push(field);rows.push(row);}
    if(!rows.length)return [];
    const headers=rows.shift();return rows.map(r=>Object.fromEntries(headers.map((h,i)=>[h,r[i]??''])));
  }
  document.getElementById('visor-importar').addEventListener('click',()=>document.getElementById('visor-archivo-csv').click());
  document.getElementById('visor-archivo-csv').addEventListener('change',async ev=>{
    let added=0,ignored=0;
    try{
      for(const file of ev.target.files){
        for(const input of parseCSV(await file.text())){
          const c=byId.get(input.caso_B);
          const vals=['x_utm18s_m','y_utm18s_m','elevacion_dem_m','pendiente_grados'].map(k=>
            input[k]==null||String(input[k]).trim()===''?NaN:Number(input[k]));
          if(!c || !vals.every(Number.isFinite) || input.crs!=='EPSG:32718' ||
             (input.rid && Number(input.rid)!==c.cfg.rid) || (input.mis && input.mis!==c.cfg.mis)){ignored++;continue;}
          const current=saved(c);
          if(current.some(p=>Math.abs(p.x_utm18s_m-vals[0])<1e-5 && Math.abs(p.y_utm18s_m-vals[1])<1e-5 &&
            Math.abs(p.elevacion_dem_m-vals[2])<1e-5 && p.tipo===(input.tipo||'B candidato') &&
            p.etiqueta===(input.etiqueta||''))){ignored++;continue;}
          current.push({...input,id_punto:current.reduce((m,p)=>Math.max(m,Number(p.id_punto)||0),0)+1,
            rid:c.cfg.rid,rio:c.cfg.rio,mis:c.cfg.mis,caso_B:c.cfg.caso_B,tipo:input.tipo||'B candidato',
            etiqueta:input.etiqueta||'',nota:input.nota||'',x_utm18s_m:vals[0],y_utm18s_m:vals[1],
            elevacion_dem_m:vals[2],pendiente_grados:vals[3],fecha_utc:input.fecha_utc||new Date().toISOString()});
          localStorage.setItem(keyOf(c.cfg),JSON.stringify({points:current,caso_B:c.cfg.caso_B}));added++;
        }
      }
      await gd.swathPuntos.recargarPuntos();counts();writeGlobal();
      message.textContent=added+' punto(s) importado(s); '+ignored+' duplicados o incompatibles omitidos.';
    }catch(e){message.textContent='No se completó la importación: '+e.message;}
    ev.target.value='';
  });

  document.getElementById('visor-color-pendiente').addEventListener('change',()=>document.getElementById('visor-aplicar').click());
  let initial;
  try{initial=byId.get(localStorage.getItem('swath_visor_caso_activo_v1'));}catch(_){}
  initial=initial||cases.find(c=>c.cfg.rid===9)||cases[0];
  width.value=initial.cfg.ancho_ew_m;length.value=initial.cfg.largo_ns_m;ev.value=initial.cfg.exageracion_vertical;
  width.max=initial.raster.max_ancho_m;length.max=initial.raster.max_largo_m;
  const first=SwathVista.ajustarVista(initial,parameters());
  rio.value=String(initial.cfg.rid);fillCases(rio.value,initial.cfg.caso_B);
  lock(true);
  await Plotly.newPlot(gd,first.fig.data,first.fig.layout,dataset.plotConfig);
  contextWrap.hidden=!first.contexto;
  if(first.contexto)await Plotly.newPlot(contextGD,first.contexto.data,first.contexto.layout,dataset.plotConfig);
  await applyTerraces();terraceInfo(first);gridInfo.textContent=first.info;
  initial={...initial,cfg:first.cfg};
  // El selector se instala una vez y recibe un contexto nuevo al cambiar de escenario.
  __INSTALAR_SELECTOR__
  active=initial;remembered.set(String(initial.cfg.rid),initial.cfg.caso_B);
  counts();message.textContent=initial.cfg.rio+' · '+initial.etiqueta+' · Ctrl+clic para marcar puntos.';
  lock(false);
})();
