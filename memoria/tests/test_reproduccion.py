"""Regression tests against final saved thesis results, not generated fixtures."""
import unittest
from pathlib import Path
import sys
import numpy as np

BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from src import analisis as a


class ResultadosFinales(unittest.TestCase):
    def test_tablas_finales(self):
        result=a.verificar()
        self.assertEqual(result['escenarios_tiempo'],27)
        self.assertEqual(result['pt_fisicos'],18)

    def test_c2_y_alternativas(self):
        d=a.calibrar_reloca()
        self.assertTrue((~d.query("frontera == 'C1'").admisible).all())
        np.testing.assert_allclose(d.query("frontera == 'C2' and mis_donante == '5e'").tau_kp_ka,[710.0539001926887],atol=1e-7)
        np.testing.assert_allclose(d.query("frontera == 'C3' and mis_donante == '5e'").tau_kp_ka,[570.8655113873691],atol=1e-7)

    def test_rahue_referencia_final(self):
        d=a.tiempos_finales().query('rid == 8')
        np.testing.assert_allclose(d.tau_G_ka,[413.82],atol=.005,rtol=0)

    def test_estado_potencial_y_ausencias(self):
        r=a.leer('catalogo_rios.csv').set_index('rid')
        self.assertEqual(r.loc[7,'estado_kp'],'sin_kp_aceptado')
        self.assertFalse(r.loc[1,'tiene_K_adoptado'])
        self.assertEqual(set(r.query("estado_kp == 'potencial'").index),{11,14,17})

    def test_identificadores_no_cambian_la_posicion(self):
        k=a.leer('catalogo_kp.csv');e=a.leer('equivalencias_identificadores_kp.csv')
        self.assertTrue(k.id_kp.is_unique)
        self.assertEqual(set(k.id_kp),set(e.id_kp))
        self.assertEqual(k.query('rid == 14').id_fuente_actual.iloc[0],'R14_ANT01')

    def test_alzamiento_sur_y_cota_absoluta(self):
        u=a.alzamiento_tramo().query("metodo == 'G' and mis == '5e' and rid in @a.SOUTH")
        np.testing.assert_allclose([u.U_KP_mm_anio.min(),u.U_KP_mm_anio.max()],[.174,.463],atol=.0005,rtol=0)
        self.assertTrue((u.cota_absoluta_sobre_tau_mm_anio>=u.U_KP_mm_anio).all())


if __name__=='__main__':unittest.main()
