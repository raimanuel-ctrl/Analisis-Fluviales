
import sys
from pathlib import Path
BASE = Path(__file__).resolve()
BASE = next(p for p in BASE.parents if (p / 'Configuracion.json').is_file())
from pathlib import Path
from datetime import datetime
import json,sys
from PIL import Image
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.pagesizes import A4,A3,landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.utils import simpleSplit,ImageReader

def crear(out):
 out=Path(out);entries=json.loads((out/'Catalogo_figuras.json').read_text(encoding='utf-8'))
 regular=Path('C:/Windows/Fonts/arial.ttf');bold=Path('C:/Windows/Fonts/arialbd.ttf')
 pdfmetrics.registerFont(TTFont('Informe',str(regular)));pdfmetrics.registerFont(TTFont('InformeBold',str(bold)))
 file=out/'Informe_figuras_20261020.pdf';canvas=Canvas(str(file),pagesize=A4,pageCompression=1)
 canvas.setTitle('Informe de figuras de geomorfología fluvial');canvas.setAuthor('Raimundo Sanchez Mosqueira')
 total=len(entries)+2
 def footer(number,W,H):
  canvas.setFillColorRGB(.35,.35,.35);canvas.setFont('Informe',8)
  canvas.drawString(32,20,'Geomorfología fluvial · θ = 0,40 · 20 octubre 2026');canvas.drawRightString(W-32,20,f'{number} / {total}');canvas.setFillColorRGB(0,0,0)
 def paragraph(text,y,size=10,width=520,x=38,leading=15):
  canvas.setFont('Informe',size)
  for line in simpleSplit(text,'Informe',size,width):canvas.drawString(x,y,line);y-=leading
  return y-8
 W,H=A4;canvas.setFont('InformeBold',20);canvas.drawString(38,H-60,'Informe de figuras fluviales')
 y=H-95;y=paragraph('Perfiles de 18 cauces, calibración de K y tiempos vigentes de los KP',y,13,width=W-76,leading=18)
 y=paragraph('20 de octubre de 2026',y,11,width=W-76)
 y=paragraph(f'El informe reúne {len(entries)} figuras con θ = 0,40, A₀ = 10⁶ m² y n = 1. Comienza con perfiles completos hasta la costa y el recorte original de planicie. Incluye los perfiles longitudinales y χ, ksn local, los dos métodos de K, controles de B y Pₜ, tiempos de los KP, comparaciones por MIS, conexiones norte–sur, cuatro cortes alternativos de Reloca y la red regional desde 0,5 km² con ksn de 500 y 1000 m. Cierra con las tablas de K promedio y tau de los KP.',y,width=W-76)
 y=paragraph('Se conservan ocho casos oficiales B/MIS y tres etiquetas de Pₜ por caso. Las etiquetas coincidentes no son mediciones independientes. Los escenarios MIS 5c y 5e permanecen separados. Las edades de highstand asignadas en la Tabla 2 publicada de Jara-Muñoz et al. (2015) son 125 ± 5 ka para 5e y 106 ± 5 ka para 5c.',y,width=W-76)
 y=paragraph('El método de propagación usa K = Δχ/(A₀^θ T). El método de pendientes usa K = U_terraza/ksn, con unidades convertidas a metros por año y ksn del ajuste OLS altura–χ entre B y Pₜ. U_terraza se obtiene de la altura del catálogo geomorfológico y el nivel marino de referencia; su cronología es la publicada.',y,width=W-76)
 y=paragraph('Los τ son tiempos de propagación condicionados por el modelo. No se interpretan como dataciones tectónicas directas. U calculado por propagacion es aparente; U del método de pendientes se impone desde la terraza. La conservación de τ(Pₜ)=T es un control algebraico, no una validación independiente.',y,width=W-76)
 y=paragraph('Las franjas litológicas representan cobertura cartográfica y no verifican la roca del lecho. Las zonas sin cartografía se muestran como Sin cobertura. Las fallas punteadas corresponden a cruces geométricos con trazas CHAF; las líneas preliminares sin precisión métrica no se convierten en cruces precisos.',y,width=W-76)
 y=paragraph('Chanco no tiene KP manual en el registro utilizado: se conserva su perfil τ, sin asignar una edad a un KP automático. Las transferencias R02/R03 desde R04 y R10–R18 desde R09 son experimentales y toman el extremo inferior del perfil receptor como tau cero. Los calibrados conservan B. R01 queda sin cálculo. El CSV manual fue enriquecido conservando sus columnas originales y un respaldo.',y,width=W-76)
 footer(1,W,H);canvas.showPage()
 canvas.setPageSize(A4);canvas.setFont('InformeBold',18);canvas.drawString(38,H-55,'Contenido y fuentes')
 y=H-95;categories={}
 for e in entries:categories.setdefault(e['categoria'],[]).append(e)
 page=3
 for cat,values in categories.items():
  canvas.setFont('InformeBold',10);name=cat.replace('_',' ');canvas.drawString(38,y,name);canvas.setFont('Informe',10);canvas.drawRightString(W-38,y,f'{len(values)} figuras · pág. {page}')
  y-=27;page+=len(values)
 y-=8;y=paragraph('Fuentes y parámetros completos: Parametros_figuras.json, Catalogo_figuras.json y las tablas de Fuentes_vigentes_figuras_20261020. Cada figura conserva un ID de categoría y río. La versión HTML permite filtrar por río y abrir los originales PNG y SVG.',y,width=W-76)
 y=paragraph('Referencias: Perron y Royden (2013), doi:10.1002/esp.3302; Royden y Perron (2013), doi:10.1002/jgrf.20031; Jara-Muñoz et al. (2015), doi:10.1016/j.quascirev.2015.01.005.',y,width=W-76)
 footer(2,W,H);canvas.showPage()
 for i,e in enumerate(entries,3):
  image=Image.open(out/e['png']);iw,ih=image.size
  dense=e['categoria'].split('_')[0] in ['01','04','06','10','12','13']
  page=A3 if dense else A4
  W,H=page if ih/iw>1.05 else landscape(page);canvas.setPageSize((W,H))
  title=simpleSplit(e['titulo'],'InformeBold',12,W-64)
  canvas.setFont('InformeBold',12);y=H-35
  for line in title:canvas.drawString(32,y,line);y-=16
  canvas.setFont('Informe',8);canvas.setFillColorRGB(.35,.35,.35);canvas.drawString(32,y,e['categoria'].replace('_',' '));canvas.setFillColorRGB(0,0,0);y-=12
  cap=simpleSplit(e['leyenda'],'Informe',9,W-64);capheight=len(cap)*12
  available_height=y-55-capheight;available_width=W-64;scale=min(available_width/iw,available_height/ih)
  width,height=iw*scale,ih*scale;x=(W-width)/2;bottom=y-height
  canvas.drawImage(ImageReader(image),x,bottom,width=width,height=height,mask='auto')
  canvas.setFont('Informe',9);ty=bottom-14
  for line in cap:canvas.drawString(32,ty,line);ty-=12
  if ty<30:raise ValueError('Leyenda fuera de página: '+e['nombre'])
  footer(i,W,H);canvas.showPage()
 canvas.save();print(json.dumps({'pdf':str(file),'paginas':total,'figuras':len(entries)},ensure_ascii=False),flush=True)
 return file

if __name__=='__main__':crear(sys.argv[1])
