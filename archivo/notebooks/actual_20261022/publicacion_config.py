"""Root discovery and a read-only preflight for the public research snapshot."""
from pathlib import Path
import json

def base():
    return Path(__file__).resolve().parent

def revisar_entradas():
    root=base()
    cfg=json.loads((root/'Configuracion.json').read_text(encoding='utf8'))
    keys=['dem','kp_manual_csv','terrazas_excel','revision_muestras_excel',
          'seleccion_excel','geologia_integrada','chaf','outlets_csv','graben_shp']
    return [{'entrada':key,'archivo':str(root/cfg[key]),'disponible':(root/cfg[key]).is_file()} for key in keys]

if __name__=='__main__':
    for item in revisar_entradas():
        print(('OK ' if item['disponible'] else 'EXTERNO ')+item['entrada']+' : '+item['archivo'])
