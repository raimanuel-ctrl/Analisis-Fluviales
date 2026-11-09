(function () {
  'use strict';
  const gd = document.getElementById('{plot_id}');
  let cfg = __CONFIG_SELECCION__;
  const keyOf = c => 'swath_puntos_v1:' + c.caso_B + ':' + c.x_B.toFixed(3) + ':' + c.y_B.toFixed(3);
  let storeKey = keyOf(cfg);
  let points = [];
  let csvHandle = null;
  const handlesByCase = new Map();
  let writeQueue = Promise.resolve();
  let surfaceHover = null;
  let gesture = null;
  let gestureSequence = 0;
  const consumed = new Set();
  let selectedTrace = null;
  const columns = ['id_punto', 'caso_B', 'rid', 'rio', 'mis', 'tipo', 'etiqueta', 'nota',
    'x_utm18s_m', 'y_utm18s_m', 'elevacion_dem_m', 'pendiente_grados', 'crs',
    'fecha_utc', 'dem_fuente', 'metodo_posicion'];

  const panel = document.createElement('section');
  panel.id = 'seleccion-swath';
  panel.setAttribute('aria-label', 'Guardar puntos del relieve');
  panel.innerHTML = `<style>
    #seleccion-swath{font:14px system-ui,sans-serif;background:#f1f6fa;border:1px solid #ccdce7;
      margin:8px 12px;padding:12px 16px;border-radius:8px;color:#19374b;line-height:1.45}
    #seleccion-swath .fila{display:flex;gap:9px;align-items:center;flex-wrap:wrap;margin:8px 0}
    #seleccion-swath label{font-size:13px;display:flex;gap:6px;align-items:center}
    #seleccion-swath input,#seleccion-swath select{font:inherit;border:1px solid #aebfcc;border-radius:4px;padding:5px 7px}
    #seleccion-swath input{width:160px} #seleccion-swath button{font:inherit;padding:7px 11px;border:1px solid #91aec4;
      background:white;color:#164b70;border-radius:5px;cursor:pointer}
    #seleccion-swath button.principal{background:#176c9e;color:white;border-color:#176c9e}
    #seleccion-swath button:disabled{opacity:.5;cursor:default}
    #seleccion-swath .nota{font-size:12px;color:#516779}#swath-status{font-weight:600}
    #seleccion-swath details{margin-top:8px}#seleccion-swath table{border-collapse:collapse;font-size:12px;width:100%;margin-top:8px}
    #seleccion-swath th,#seleccion-swath td{padding:5px 7px;border-bottom:1px solid #d3dfe7;text-align:left}
    #seleccion-swath .tabla{overflow:auto;max-height:240px;background:#fff}
  </style>
  <strong>Ctrl + clic sobre el relieve 3D: marcar un punto</strong>
  <span class="nota"> · Giro normal sin Ctrl · EPSG:32718</span>
  <div id="swath-contexto" class="nota"></div>
  <div class="fila">
    <label>Tipo <select id="swath-tipo"><option>B candidato</option><option>Borde de terraza</option>
      <option>Paleoacantilado</option><option>Plataforma</option><option>Otro</option></select></label>
    <label>Etiqueta <input id="swath-etiqueta" placeholder="Opcional" maxlength="120"></label>
    <label>Nota <input id="swath-nota" placeholder="Observación" maxlength="500"></label>
    <button id="swath-vincular" class="principal">Elegir CSV y guardar automáticamente</button>
    <button id="swath-descargar">Descargar CSV</button>
    <button id="swath-deshacer">Deshacer último punto</button>
  </div>
  <div id="swath-status" role="status" aria-live="polite">Sin puntos marcados.</div>
  <div class="nota">Para escribir cada punto en el mismo archivo: elige el CSV una vez y conserva esta página abierta.
    Sin vincular un archivo, los puntos se conservan en este navegador y puedes usar «Descargar CSV».
    El clic toma un nodo de la malla del DEM; no implica precisión subpíxel. B sigue siendo un candidato.</div>
  <details id="swath-details"><summary>Ver puntos marcados (<span id="swath-count">0</span>)</summary>
    <div class="tabla"><table><thead><tr><th>ID</th><th>Tipo</th><th>Etiqueta</th><th>Este UTM (m)</th>
      <th>Norte UTM (m)</th><th>Z DEM (m)</th><th>Pendiente (°)</th></tr></thead><tbody id="swath-rows"></tbody></table></div>
  </details>`;
  gd.parentNode.insertBefore(panel, gd);
  const el = id => panel.querySelector('#' + id);
  const status = message => { el('swath-status').textContent = message; };
  const escapeCsv = v => '"' + String(v == null ? '' : v).replace(/"/g, '""') + '"';
  function csvText() {
    return '\ufeff' + columns.join(',') + '\r\n' + points.map(p =>
      columns.map(c => escapeCsv(p[c])).join(',')).join('\r\n') + (points.length ? '\r\n' : '');
  }
  function persist() {
    try { localStorage.setItem(storeKey, JSON.stringify({points: points, caso_B: cfg.caso_B})); }
    catch (_) { status('Punto marcado; este navegador no permite respaldo local. Descarga el CSV.'); }
    if (typeof window.dispatchEvent === 'function' && typeof CustomEvent === 'function')
      window.dispatchEvent(new CustomEvent('swath-puntos-cambiaron'));
  }
  function drawPoints() {
    const rows = el('swath-rows');
    rows.replaceChildren();
    points.forEach(p => {
      const tr = document.createElement('tr');
      [p.id_punto, p.tipo, p.etiqueta, Number(p.x_utm18s_m).toFixed(3), Number(p.y_utm18s_m).toFixed(3),
        Number(p.elevacion_dem_m).toFixed(2), Number(p.pendiente_grados).toFixed(2)].forEach(value => {
          const td = document.createElement('td'); td.textContent = value; tr.appendChild(td);
        });
      rows.appendChild(tr);
    });
    el('swath-count').textContent = String(points.length);
    el('swath-deshacer').disabled = !points.length;
    // Fuera del recorte se ocultan los símbolos; las filas y CSV conservan todos los puntos.
    const inside=p=>Math.abs(p.x_utm18s_m-cfg.x_B)<=(cfg.ancho_ew_m??Infinity)/2&&Math.abs(p.y_utm18s_m-cfg.y_B)<=(cfg.largo_ns_m??Infinity)/2;
    const x = points.map(p => inside(p)?p.x_utm18s_m - cfg.x_B:null);
    const y = points.map(p => inside(p)?p.y_utm18s_m - cfg.y_B:null);
    const z = points.map(p => inside(p)?p.elevacion_dem_m:null);
    const text = points.map(p => p.etiqueta || ('P' + p.id_punto));
    if (selectedTrace === null) {
      selectedTrace = gd.data.length;
      return Plotly.addTraces(gd, {type:'scatter3d', scene:'scene', x:x, y:y, z:z,
        mode:'markers+text', text:text, textposition:'top center', name:'Puntos marcados',
        marker:{size:6,color:'#fff300',line:{color:'#1d2227',width:2},symbol:'diamond'},
        hovertemplate:'%{text}<br>ΔE %{x:.1f} m · ΔN %{y:.1f} m<br>Z %{z:.2f} m<extra>Punto candidato</extra>'});
    } else {
      return Plotly.restyle(gd, {x:[x], y:[y], z:[z], text:[text]}, [selectedTrace]);
    }
  }
  function saveToFile() {
    if (!csvHandle) return;
    // Capturar archivo y filas antes de un posible cambio de río.
    const handle = csvHandle, text = csvText(), key = storeKey, count = points.length;
    writeQueue = writeQueue.then(async () => {
      const writable = await handle.createWritable();
      try { await writable.write(text); await writable.close(); }
      catch (e) { try { await writable.abort(); } catch (_) {} throw e; }
      if (storeKey === key) status(count + ' punto(s). CSV actualizado: ' + handle.name);
    }).catch(e => {
      handlesByCase.delete(key);
      if (storeKey === key) {
        csvHandle = null;
        status('No se pudo actualizar el CSV. Los puntos siguen en la página: usa «Descargar CSV».');
      }
      console.warn('Escritura CSV:', e.message);
    });
  }
  function pointFromEvent(ev) {
    const list = ev && ev.points;
    if (!list) return null;
    const p = list.find(p => p.curveNumber === 0 && p.data && p.data.type === 'surface');
    if (!p) return null;
    const ex = Number(p.x), ny = Number(p.y);
    if (![ex,ny].every(Number.isFinite)) return null;
    const ix = Math.round((ex-cfg.east_m[0])/(cfg.east_m[1]-cfg.east_m[0]));
    const iy = Math.round((ny-cfg.north_m[0])/(cfg.north_m[1]-cfg.north_m[0]));
    if (ix < 0 || iy < 0 || ix >= cfg.east_m.length || iy >= cfg.north_m.length) return null;
    const E = cfg.x_B + cfg.east_m[ix], N = cfg.y_B + cfg.north_m[iy];
    const Z = cfg.z_m[iy][ix], slope = cfg.pendiente_grados[iy][ix];
    if (![E, N, Z, slope].every(Number.isFinite)) return null;
    return {E:E, N:N, Z:Z, slope:slope};
  }
  function addPoint(p, id) {
    if (!p || consumed.has(id)) return;
    consumed.add(id);
    const maxid = points.reduce((m, v) => Math.max(m, Number(v.id_punto) || 0), 0);
    points.push({id_punto:maxid+1, caso_B:cfg.caso_B, rid:cfg.rid, rio:cfg.rio || '', mis:cfg.mis,
      tipo:el('swath-tipo').value, etiqueta:el('swath-etiqueta').value.trim(), nota:el('swath-nota').value.trim(),
      x_utm18s_m:p.E, y_utm18s_m:p.N, elevacion_dem_m:p.Z, pendiente_grados:p.slope,
      crs:'EPSG:32718', fecha_utc:new Date().toISOString(), dem_fuente:cfg.dem,
      metodo_posicion:'Ctrl+clic sobre nodo de malla 3D; candidato sin validar'});
    persist(); drawPoints();
    status(points.length + ' punto(s). ' + (csvHandle ? 'Guardando CSV…' : 'Respaldo local; vincula o descarga el CSV.'));
    saveToFile();
  }
  gd.on('plotly_hover', ev => { surfaceHover = pointFromEvent(ev); });
  gd.on('plotly_unhover', () => { surfaceHover = null; });
  gd.addEventListener('pointerdown', ev => {
    if (ev.button !== 0) return;
    gesture = {id:++gestureSequence, key:storeKey, ctrl:ev.ctrlKey || ev.metaKey, x:ev.clientX, y:ev.clientY,
      started:performance.now(), point:surfaceHover, valid:false};
  }, true);
  gd.addEventListener('pointerup', ev => {
    if (!gesture) return;
    const g = gesture;
    g.valid = g.ctrl && Math.hypot(ev.clientX-g.x,ev.clientY-g.y) < 7 && performance.now()-g.started < 1500;
    if (!g.valid) return;
    g.point = surfaceHover || g.point;
    // plotly_click aporta la selección de superficie; hover es respaldo si el motor 3D no emite click.
    setTimeout(() => { if (g.key === storeKey && !consumed.has(g.id)) addPoint(g.point, g.id); }, 200);
  }, true);
  gd.on('plotly_click', ev => {
    if (gesture && gesture.valid) addPoint(pointFromEvent(ev), gesture.id);
  });
  el('swath-vincular').addEventListener('click', async () => {
    if (!window.showSaveFilePicker) {
      status('Este navegador no permite escritura directa. Usa Ctrl+clic y luego «Descargar CSV».');
      return;
    }
    try {
      csvHandle = await window.showSaveFilePicker({suggestedName:'Puntos_' + cfg.caso_B + '.csv',
        types:[{description:'Tabla CSV de puntos',accept:{'text/csv':['.csv']}}]});
      handlesByCase.set(storeKey, csvHandle);
      status('CSV vinculado: ' + csvHandle.name + '. Ctrl+clic guardará cada punto.');
      saveToFile();
    } catch (e) { if (e.name !== 'AbortError') status('No se pudo vincular el archivo; usa «Descargar CSV».'); }
  });
  el('swath-descargar').addEventListener('click', () => {
    const url = URL.createObjectURL(new Blob([csvText()], {type:'text/csv;charset=utf-8;'}));
    const a = document.createElement('a'); a.href = url; a.download = 'Puntos_' + cfg.caso_B + '.csv';
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 15000);
    status(points.length + ' punto(s) exportado(s) a CSV.');
  });
  el('swath-deshacer').addEventListener('click', () => {
    if (!points.length) return;
    points.pop(); persist(); drawPoints();
    status(points.length + ' punto(s). Último punto retirado.'); saveToFile();
  });
  function restore() {
    points = [];
    try {
      const saved = JSON.parse(localStorage.getItem(storeKey) || 'null');
      if (saved && Array.isArray(saved.points)) points = saved.points.filter(p =>
        p.caso_B === cfg.caso_B && [p.x_utm18s_m,p.y_utm18s_m,p.elevacion_dem_m,p.pendiente_grados].every(Number.isFinite))
          .map(p => ({...p, rio:cfg.rio || p.rio || ''}));
    } catch (_) {}
    el('swath-contexto').textContent = (cfg.rio || ('R' + cfg.rid)) + ' · MIS ' + cfg.mis + ' · ' + cfg.caso_B;
  }
  gd.swathPuntos = {
    async prepararCambio() {
      if (gesture && gesture.valid && gesture.key === storeKey && !consumed.has(gesture.id))
        addPoint(gesture.point, gesture.id);
      await writeQueue;
      gesture = null; surfaceHover = null;
    },
    async cambiarCaso(next) {
      const different=keyOf(next)!==storeKey;
      cfg = next; storeKey = keyOf(cfg); csvHandle = handlesByCase.get(storeKey) || null;
      selectedTrace = null; gesture = null; surfaceHover = null;
      if(different){el('swath-etiqueta').value = ''; el('swath-nota').value = '';}
      restore(); await drawPoints();
      status(points.length + ' punto(s) en este escenario. ' + (csvHandle ? 'CSV vinculado.' : 'Vincula o descarga su CSV.'));
    },
    async recargarPuntos() { restore(); await drawPoints(); saveToFile(); },
    puntos() { return points.map(p => ({...p})); }
  };
  restore();
  drawPoints();
  status(points.length ? points.length + ' punto(s) recuperado(s). Vuelve a vincular el CSV para escritura automática.' :
    'Listo: Ctrl+clic en el relieve. Puedes elegir primero el archivo CSV.');
})();
