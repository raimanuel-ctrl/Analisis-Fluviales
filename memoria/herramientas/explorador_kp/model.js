(function (root) {
'use strict';
const FACTOR=Math.pow(1e6,.4), NOW=()=>new Date().toISOString();
function finite(x){return typeof x==='number'&&Number.isFinite(x)}
function uuid(){return globalThis.crypto?.randomUUID?.()||'kp-'+Date.now().toString(36)+'-'+Math.random().toString(36).slice(2)}
function tau(chi,K,origin=0){if(!(finite(K)&&K>0))throw Error('K debe ser positivo');return (chi-origin)/(K*FACTOR*1000)}
function lowerBound(a,v,get=x=>x){let l=0,r=a.length;while(l<r){let m=(l+r)>>1;if(get(a[m])<v)l=m+1;else r=m}return l}
function point(arr,i,pid,tid=null){return {s:arr[0][i],chi:arr[1][i],z:arr[2][i],area:arr[3][i],x:arr[4][i],y:arr[5][i],pid,tid,index:i}}
function interp(points,value,key='chi'){
 if(!points.length||value<points[0][key]-1e-7||value>points.at(-1)[key]+1e-7)return null;
 let i=lowerBound(points,value,p=>p[key]);if(i===0)return {...points[0]};if(i>=points.length)return {...points.at(-1)};
 const a=points[i-1],b=points[i],f=(value-a[key])/(b[key]-a[key]||1),q={...((f<.5)?a:b),interpolated:true};for(const k of ['s','chi','z','area','x','y'])q[k]=a[k]+f*(b[k]-a[k]);return q;
}
function ols(points){
 if(points.length<5)return null;let mx=0,my=0;for(const p of points){mx+=p.chi;my+=p.z}mx/=points.length;my/=points.length;let xx=0,yy=0,xy=0;
 for(const p of points){let x=p.chi-mx,y=p.z-my;xx+=x*x;yy+=y*y;xy+=x*y}if(xx<1e-12)return null;
 let b=xy/xx,a=my-b*mx,r2=yy>1e-18?Math.min(1,Math.max(0,xy*xy/(xx*yy))):null,ss=0;for(const p of points)ss+=(p.z-(a+b*p.chi))**2;
 return {slope:b,intercept:a,R2:r2,RMSE:Math.sqrt(ss/points.length),n:points.length,span:points.at(-1).s-points[0].s,chiMin:points[0].chi,chiMax:points.at(-1).chi,sMin:points[0].s,sMax:points.at(-1).s};
}
class Network{
 constructor(desc){this.desc=desc;this.pid=desc.id;this.main=Array.from({length:desc.main[0].length},(_,i)=>point(desc.main,i,desc.id));this.tributaries=(desc.tributaries||[]).map(t=>({id:t.id,points:Array.from({length:t.data[0].length},(_,i)=>point(t.data,i,desc.id,t.id)),sConfluence:t.s_confluence,chiConfluence:t.chi_confluence}));this.trById=new Map(this.tributaries.map(t=>[t.id,t]));this.nodes=new Map();this.routes=new Map();
  const key=p=>p.x.toFixed(3)+','+p.y.toFixed(3);this.key=key;
  // A common spatial node keeps the main-channel values when available.
  for(const t of this.tributaries)for(const p of t.points)if(!this.nodes.has(key(p)))this.nodes.set(key(p),{p,down:null,up:[]});
  for(const p of this.main)this.nodes.set(key(p),{p,down:null,up:[]});
  const connect=pts=>{for(let i=1;i<pts.length;i++){const a=this.nodes.get(key(pts[i])),b=this.nodes.get(key(pts[i-1]));if(a===b||a.p.s<=b.p.s+1e-7)continue;if(!a.down||pts[i].tid===null)a.down=b;}};
  for(const t of this.tributaries)connect(t.points);connect(this.main);
  for(const n of this.nodes.values())if(n.down)n.down.up.push(n);
 }
 data(tid){return tid?this.trById.get(tid)?.points:this.main}
 exact(tid,i){const a=this.data(tid);if(!a||i<0||i>=a.length)throw Error('Punto fuera del perfil');return {...a[i]}}
 route(center){
  if(!center.tid)return this.main;
  const node=this.nodes.get(this.key(center));if(!node)return [center];let lower=[],cur=node,seen=new Set();
  while(cur&&!seen.has(cur)){seen.add(cur);lower.push(cur.p);cur=cur.down}
  lower.reverse();let upper=[];cur=node;seen=new Set([node]);
  while(cur.up.length){let next=cur.up.reduce((a,b)=>(b.p.area>a.p.area?b:a));if(seen.has(next))break;seen.add(next);upper.push(next.p);cur=next}
  return [...lower,...upper];
 }
 excludedByCut(p){const c=this.desc.regression_cut;if(!c||p.s>=c.s_min_m-1e-7)return false;return p.tid===null||this.nodes.get(this.key(p))?.p.tid===null}
 regression(center,W='full',gap=100){
  const route=this.route(center),retained=route.filter(p=>!this.excludedByCut(p)),full=W==='full',lo=retained.filter(p=>(full||p.s>=center.s-W-1e-7)&&p.s<center.s-gap),hi=retained.filter(p=>p.s>center.s+gap&&(full||p.s<=center.s+W+1e-7)),lower=ols(lo),upper=ols(hi),reviewed=full?retained:retained.filter(p=>p.s>=center.s-W&&p.s<=center.s+W),result={center:{...center},W:full?null:W,scope:full?'perfil_completo':'ventana_fija',exclusionAluvial:this.desc.regression_cut||null,nExcluded:route.length-retained.length,gap,lower,upper,route,reviewedPoints:reviewed,lowerPoints:lo,upperPoints:hi,cross:null,projected:null,reason:''};
  if(this.excludedByCut(center)){result.lower=null;result.upper=null;result.reason='Punto dentro del tramo bajo excluido por el corte documentado';return result}
  if(!lower||!upper){result.reason='Faltan nodos para dos ajustes (mínimo 5 por lado)';return result}
  const den=lower.slope-upper.slope;if(Math.abs(den)<1e-12){result.reason='Las rectas son paralelas o casi paralelas';return result}
  const chi=(upper.intercept-lower.intercept)/den,z=lower.intercept+lower.slope*chi;result.cross={chi,z};result.projected=interp(reviewed,chi);
  if(!result.projected)result.reason='El cruce queda fuera del tramo revisado';return result;
 }
}
function validateRecord(r,data){
 if(!data.profiles.some(p=>p.id===r.id_perfil)||!data.catalog.some(k=>k.id===r.id_K))throw Error('Una medición usa un perfil o un K que no está en esta copia');
 if(!['x_utm_m','y_utm_m','altura_kp_m','chi_m','K_m02_anio','tau_ka'].every(k=>finite(r[k])))throw Error('Hay una medición con valores no válidos');
 const k=data.catalog.find(k=>k.id===r.id_K),o=data.profiles.find(p=>p.id===r.id_perfil).origins.find(o=>o.id===r.id_origen);if(!o)throw Error('Origen temporal no disponible');
 if(Math.abs(k.K_m02_anio-r.K_m02_anio)>Math.max(1e-18,Math.abs(k.K_m02_anio)*1e-10))throw Error('El K de la medición no coincide con su catálogo');
 if(Math.abs(tau(r.chi_m,r.K_m02_anio,o.chi_m)-r.tau_ka)>Math.max(1e-7,Math.abs(r.tau_ka)*1e-9))throw Error('No se reproduce el tiempo de una medición');return true;
}
function normalizeState(raw,data){
 let s=raw.state||raw;if(raw.dataset_id&&raw.dataset_id!==data.id)throw Error('Este proyecto corresponde a otra copia de datos');
 if(raw.parametros&&(raw.parametros.theta!==.4||raw.parametros.A0_m2!==1e6||raw.parametros.n!==1))throw Error('Los parámetros del proyecto no corresponden a θ=0,40, A0=1 km² y n=1');
 s=JSON.parse(JSON.stringify(s));s.mediciones=s.mediciones||[];s.contadores=s.contadores||{};s.historial=s.historial||[];s.notas_proyecto=s.notas_proyecto||'';s.vista=s.vista||{};for(const r of s.mediciones)validateRecord(r,data);return s;
}
function addRecord(state,data,network,p,k,origin,options={}){
 const desc=network.desc;
 if(p.tid){const tr=network.trById.get(p.tid);if(tr&&origin.s_m>tr.sConfluence+1e-7)throw Error('El origen está aguas arriba de la confluencia de este tributario');}
 const elapsed=tau(p.chi,k.K_m02_anio,origin.chi_m);if(elapsed< -1e-9)throw Error('El punto está aguas abajo del origen de τ');
 const existing=state.mediciones.find(r=>r.id_perfil===desc.id&&r.id_K===k.id&&r.id_origen===origin.id&&Math.hypot(r.x_utm_m-p.x,r.y_utm_m-p.y)<1e-5);if(existing)return {record:existing,created:false};
 const same=state.mediciones.find(r=>r.id_rio===desc.rid&&Math.hypot(r.x_utm_m-p.x,r.y_utm_m-p.y)<1e-5);let id=same?.id_kp;if(!id){let count=(state.contadores[desc.rid]||0)+1;while(state.mediciones.some(r=>r.id_kp===`R${String(desc.rid).padStart(2,'0')}_KP${String(count).padStart(3,'0')}`))count++;state.contadores[desc.rid]=count;id=`R${String(desc.rid).padStart(2,'0')}_KP${String(count).padStart(3,'0')}`;}
 const r={id_medicion:uuid(),id_rio:desc.rid,rio:desc.rio,id_kp:id,origen_K:k.label,id_K:k.id,K_m02_anio:k.K_m02_anio,metodo_K:k.metodo,rio_donante:k.rid_donante,MIS_K:k.mis,geometria_donante:k.geometria_donante,altura_kp_m:p.z,tau_ka:Math.max(0,elapsed),distancia_longitudinal_m:p.s,chi_m:p.chi,x_utm_m:p.x,y_utm_m:p.y,area_m2:p.area,id_perfil:desc.id,geometria_receptor:desc.geometry,tipo_perfil_clic:options.panel||'HTML',indice_nodo:p.index,metodo_seleccion:p.interpolated?'interseccion_regresion_chi_proyectada_en_perfil':'nodo mas cercano en pantalla',id_origen:origin.id,origen_tau:origin.label,s_origen_m:origin.s_m,chi_origen_m:origin.chi_m,delta_chi_m:p.chi-origin.chi_m,theta:.4,A0_m2:1e6,n:1,crs:'EPSG:32718',perfil_sha256:desc.sha256,fecha_utc:NOW(),nota:options.note||'',estado:'exploratorio',tipo_cauce:p.tid?'tributario':'principal',id_tributario:p.tid||null};
 if(p.tid){let tr=network.trById.get(p.tid);r.s_confluencia_principal_m=tr?.sConfluence;r.chi_confluencia_principal_m=tr?.chiConfluence;}
 if(options.regression){const a=options.regression;r.regresion={corte_documentado:a.exclusionAluvial,n_nodos_excluidos_por_corte:a.nExcluded,alcance:a.scope,limite_inferior_s_m:a.reviewedPoints[0]?.s,limite_superior_s_m:a.reviewedPoints.at(-1)?.s,ventana_por_lado_m:a.W,exclusion_por_lado_m:a.gap,centro_s_m:a.center.s,centro_chi_m:a.center.chi,inferior:a.lower,superior:a.upper,chi_interseccion_m:a.cross?.chi??null,z_rectas_m:a.cross?.z??null,z_perfil_m:a.projected?.z??null,solo_exploratoria:true};}
 state.mediciones.push(r);state.historial.push({fecha_utc:NOW(),accion:'crear',id_medicion:r.id_medicion});state.actualizado_utc=NOW();return {record:r,created:true};
}
function csv(records){const base=['id_rio','rio','id_kp','origen_K','id_K','altura_kp_m','tau_ka','distancia_longitudinal_m','chi_m','x_utm_m','y_utm_m','tipo_cauce','id_tributario','K_m02_anio','MIS_K','id_origen','chi_origen_m','theta','A0_m2','n','metodo_seleccion','fecha_utc','nota'];const keys=[...new Set([...base,...records.flatMap(r=>Object.keys(r))])];const esc=v=>{if(v===null||v===undefined)return '';let s=typeof v==='object'?JSON.stringify(v):String(v);return /[",\r\n]/.test(s)?'"'+s.replaceAll('"','""')+'"':s};return '\ufeff'+[keys.join(','),...records.map(r=>keys.map(k=>esc(r[k])).join(','))].join('\r\n')+'\r\n';}
const enc=new TextEncoder();
function shp(records){
 let n=records.length,buf=new ArrayBuffer(100+44*n),v=new DataView(buf),ix=new ArrayBuffer(100+8*n),iv=new DataView(ix);let xs=records.map(r=>r.x_utm_m),ys=records.map(r=>r.y_utm_m),zs=records.map(r=>r.altura_kp_m),range=a=>a.length?[Math.min(...a),Math.max(...a)]:[0,0],xr=range(xs),yr=range(ys),zr=range(zs);
 for(const [b,len] of [[v,buf.byteLength],[iv,ix.byteLength]]){b.setInt32(0,9994);b.setInt32(24,len/2);b.setInt32(28,1000,true);b.setInt32(32,11,true);[xr[0],yr[0],xr[1],yr[1],zr[0],zr[1],-1e39,-1e39].forEach((x,i)=>b.setFloat64(36+i*8,x,true));}
 records.forEach((r,i)=>{let o=100+i*44;v.setInt32(o,i+1);v.setInt32(o+4,18);v.setInt32(o+8,11,true);[r.x_utm_m,r.y_utm_m,r.altura_kp_m,-1e39].forEach((x,j)=>v.setFloat64(o+12+j*8,x,true));iv.setInt32(100+i*8,o/2);iv.setInt32(104+i*8,18);});
 const fields=[['med_id','C',36,0,'id_medicion'],['rio_id','N',5,0,'id_rio'],['kp_id','C',30,0,'id_kp'],['k_id','C',100,0,'id_K'],['k_origen','C',180,0,'origen_K'],['mis_k','C',5,0,'MIS_K'],['k_valor','N',24,16,'K_m02_anio'],['z_m','N',20,7,'altura_kp_m'],['tau_ka','N',22,9,'tau_ka'],['s_m','N',22,7,'distancia_longitudinal_m'],['chi_m','N',22,7,'chi_m'],['x_utm','N',20,6,'x_utm_m'],['y_utm','N',20,6,'y_utm_m'],['perfil','C',60,0,'id_perfil'],['origen','C',60,0,'id_origen'],['chi_0','N',22,7,'chi_origen_m'],['theta','N',5,2,'theta'],['a0_m2','N',16,1,'A0_m2'],['n_exp','N',5,2,'n'],['cauce','C',12,0,'tipo_cauce'],['trib_id','C',80,0,'id_tributario'],['nota','C',254,0,'nota']];
 const h=32+32*fields.length+1,rl=1+fields.reduce((s,f)=>s+f[2],0),db=new Uint8Array(h+rl*n+1),dv=new DataView(db.buffer),d=new Date();db.fill(32);db[0]=3;db[1]=d.getFullYear()-1900;db[2]=d.getMonth()+1;db[3]=d.getDate();dv.setUint32(4,n,true);dv.setUint16(8,h,true);dv.setUint16(10,rl,true);for(let i=12;i<32;i++)db[i]=0;
 fields.forEach((f,i)=>{let o=32+i*32;for(let j=0;j<32;j++)db[o+j]=0;db.set(enc.encode(f[0]).slice(0,10),o);db[o+11]=f[1].charCodeAt(0);db[o+16]=f[2];db[o+17]=f[3];});db[h-1]=13;
 records.forEach((r,i)=>{let o=h+i*rl;db[o++]=32;for(const f of fields){let val=r[f[4]],text='';if(f[1]==='N'){if(finite(val))text=val.toFixed(f[3]);if(text.length>f[2])throw Error('Un valor excede la precision del SHP: '+f[0]);text=text.padStart(f[2],' ');db.set(enc.encode(text),o);}else{let bytes=[];for(const c of String(val??'')){let a=[...enc.encode(c)];if(bytes.length+a.length>f[2])break;bytes.push(...a)}db.set(bytes,o);}o+=f[2];}});db[db.length-1]=26;
 const prj='PROJCS["WGS 84 / UTM zone 18S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",-75],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1],AUTHORITY["EPSG","32718"]]';return {shp:new Uint8Array(buf),shx:new Uint8Array(ix),dbf:db,prj:enc.encode(prj),cpg:enc.encode('UTF-8')};
}
const crcTable=Array.from({length:256},(_,n)=>{let c=n;for(let j=0;j<8;j++)c=(c&1)?0xedb88320^(c>>>1):c>>>1;return c>>>0});
function crc32(b){let c=0xffffffff;for(let i=0;i<b.length;i++)c=crcTable[(c^b[i])&255]^(c>>>8);return (c^0xffffffff)>>>0}
function zip(files){const locals=[],centrals=[];let offset=0;const today=new Date(),dosDate=((today.getFullYear()-1980)<<9)|((today.getMonth()+1)<<5)|today.getDate(),dosTime=(today.getHours()<<11)|(today.getMinutes()<<5)|(today.getSeconds()>>1);for(const [name,raw] of files){let bytes=typeof raw==='string'?enc.encode(raw):raw,fn=enc.encode(name),crc=crc32(bytes),h=new Uint8Array(30+fn.length),v=new DataView(h.buffer);v.setUint32(0,0x04034b50,true);v.setUint16(4,20,true);v.setUint16(6,0x800,true);v.setUint16(10,dosTime,true);v.setUint16(12,dosDate,true);v.setUint32(14,crc,true);v.setUint32(18,bytes.length,true);v.setUint32(22,bytes.length,true);v.setUint16(26,fn.length,true);h.set(fn,30);locals.push(h,bytes);let c=new Uint8Array(46+fn.length),cv=new DataView(c.buffer);cv.setUint32(0,0x02014b50,true);cv.setUint16(4,20,true);cv.setUint16(6,20,true);cv.setUint16(8,0x800,true);cv.setUint16(12,dosTime,true);cv.setUint16(14,dosDate,true);cv.setUint32(16,crc,true);cv.setUint32(20,bytes.length,true);cv.setUint32(24,bytes.length,true);cv.setUint16(28,fn.length,true);cv.setUint32(42,offset,true);c.set(fn,46);centrals.push(c);offset+=h.length+bytes.length;}let size=centrals.reduce((a,b)=>a+b.length,0),end=new Uint8Array(22),ev=new DataView(end.buffer);ev.setUint32(0,0x06054b50,true);ev.setUint16(8,files.length,true);ev.setUint16(10,files.length,true);ev.setUint32(12,size,true);ev.setUint32(16,offset,true);let all=[...locals,...centrals,end],out=new Uint8Array(all.reduce((a,b)=>a+b.length,0)),pos=0;for(const b of all){out.set(b,pos);pos+=b.length}return out;}
const API={FACTOR,NOW,finite,uuid,tau,lowerBound,point,interp,ols,Network,normalizeState,validateRecord,addRecord,csv,shp,zip,crc32};if(typeof module!=='undefined'&&module.exports)module.exports=API;root.KPModel=API;
})(typeof globalThis!=='undefined'?globalThis:this);
