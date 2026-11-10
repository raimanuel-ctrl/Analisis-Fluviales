"""Compare the plotted fit inputs with independent OLS reconstruction."""
from pathlib import Path
import json,sys,tempfile,unittest
import numpy as np
import pandas as pd
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from src import analisis as a
from src.graficos_informe import Graficos,MAPS


class GraficosInforme(unittest.TestCase):
    def test_inventario_y_mapas(self):
        d=json.loads((a.DATA/'graficos/inventario_figuras.json').read_text(encoding='utf-8'))
        self.assertEqual({r['figura'] for r in d},set(range(1,57)))
        self.assertEqual({r['figura'] for r in d if r['tipo']=='mapa_excluido'},MAPS)
        self.assertEqual(len(MAPS),8)

    def test_68_regresiones_sobre_los_nodos_originales(self):
        fits=pd.read_csv(a.DATA/'graficos/regresiones_completo_3km.csv',float_precision='round_trip');n=0
        for _,r in fits.iterrows():
            if r.rid==7:continue
            p=a.perfil(r.rid)
            for scope in ['completo','3km']:
                for side in ['inferior','superior']:
                    q=p[p.distancia_m.between(r[f'{scope}_{side}_s_min_m']-1e-6,r[f'{scope}_{side}_s_max_m']+1e-6)]
                    m,b=np.polyfit(q.chi_m,q.elevacion_m,1)
                    self.assertEqual(len(q),r[f'{scope}_{side}_n'])
                    np.testing.assert_allclose([m,b],[r[f'{scope}_{side}_b'],r[f'{scope}_{side}_a']],rtol=1e-10,atol=1e-10);n+=1
        self.assertEqual(n,68)

    def test_216_ajustes_de_incision(self):
        n=0
        for r in a.incision_ventanas().itertuples():
            p=a.perfil(r.rid);q=p[p.distancia_m.between(r.s_min_m-1e-6,r.s_max_m+1e-6)]
            m,b=np.polyfit(q.chi_m,q.elevacion_m,1)
            self.assertEqual(len(q),r.n_nodos)
            np.testing.assert_allclose(m*a.FACTOR,r.ksn_m08,rtol=1e-10,atol=1e-10);n+=1
        self.assertEqual(n,216)

    def test_entradas_externas_y_panel_mixto(self):
        temporary_root=(BASE/'resultados').resolve()
        assert temporary_root.is_relative_to(BASE.resolve())
        temporary_root.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='prueba_graficos_',dir=temporary_root) as root:
            assert Path(root).resolve().is_relative_to(temporary_root)
            g=Graficos(Path(root)/'salidas',Path(root)/'sin_entradas',formats=('svg',))
            d=g.generar([1,2,14,16,29,56]);status={r['figura']:r['estado'] for r in d}
            self.assertEqual(status[1],'mapa_excluido')
            self.assertEqual({n for n,s in status.items() if s=='requiere_entrada_local'},{2,16,29})
            self.assertEqual(status[14],'generada');self.assertEqual(status[56],'generada')
            self.assertFalse((g.out/'Figura_01.svg').exists())
            self.assertTrue((g.out/'Figura_14_solo_perfil.svg').is_file())


if __name__=='__main__':unittest.main()
