"""
Tests de saber_brr.

* Sintéticos: siempre corren, no dependen de los datos.
* Con datos reales: se omiten si no se define la variable de entorno
  SABER_DATA (ruta al .txt o .parquet de Saber 3579 2022).
  Valores esperados: sección 5.4 de la guía metodológica del ICFES.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import saber_brr as sb  # noqa: E402


# ------------------------------------------------------------- datos sintéticos
def hacer_datos(n=400, n_pv=3, n_rep=8, seed=0, w_rep_igual=False):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    d = pd.DataFrame({"peso": rng.uniform(1, 30, n), "x": x})
    for i in range(1, n_pv + 1):
        d[f"pv{i}"] = 400 + 25 * x + rng.normal(scale=60, size=n)
    for r in range(1, n_rep + 1):
        d[f"peso_r{r}"] = (d["peso"] if w_rep_igual
                           else d["peso"] * rng.choice([0.5, 1.5], size=n))
    return d


KW = dict(n_pv=3, n_rep=8)


# ------------------------------------------------------------------ sintéticos
def test_media_ponderada_pesos_iguales():
    y = np.array([1.0, 2.0, 6.0])
    assert sb.media_ponderada(y, np.ones(3)) == pytest.approx(3.0)


def test_media_ponderada_con_pesos():
    assert sb.media_ponderada([10, 20], [1, 3]) == pytest.approx(17.5)


def test_beta_wls_coincide_con_lstsq():
    d = hacer_datos()
    X = np.column_stack([np.ones(len(d)), d["x"]])
    y, w = d["pv1"].to_numpy(), d["peso"].to_numpy()
    esperado = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
    assert sb.beta_wls(X, y, w) == pytest.approx(esperado, abs=1e-9)


def test_beta_wls_coincide_con_statsmodels():
    sm = pytest.importorskip("statsmodels.api")
    d = hacer_datos()
    X = pd.DataFrame({"const": 1.0, "x": d["x"]})
    y, w = d["pv1"].to_numpy(), d["peso"].to_numpy()
    mod = sm.WLS(y, X, weights=w).fit()
    assert sb.beta_wls(X, y, w) == pytest.approx(mod.params.values, abs=1e-9)


def test_media_se_cero_sin_variacion():
    """Réplicas = peso base y PV idénticos -> error estándar 0."""
    d = hacer_datos(w_rep_igual=True)
    for i in (2, 3):
        d[f"pv{i}"] = d["pv1"]
    est, se = sb.media_pv_brr(d, "pv", **KW)
    assert est == pytest.approx(sb.media_ponderada(d["pv1"], d["peso"]))
    assert se == pytest.approx(0.0, abs=1e-9)


def test_media_solo_varianza_de_imputacion():
    """Con réplicas iguales, se = sqrt((1 + 1/M) * var entre PV)."""
    d = hacer_datos(w_rep_igual=True)
    medias = np.array([sb.media_ponderada(d[f"pv{i}"], d["peso"]) for i in (1, 2, 3)])
    esperado = np.sqrt((1 + 1 / 3) * medias.var(ddof=1))
    _, se = sb.media_pv_brr(d, "pv", **KW)
    assert se == pytest.approx(esperado)


def test_media_formula_brr_manual():
    """Un PV repetido: se = sqrt(sum((theta_r - theta)^2) / (G (1-k)^2))."""
    d = hacer_datos()
    for i in (2, 3):
        d[f"pv{i}"] = d["pv1"]
    theta = sb.media_ponderada(d["pv1"], d["peso"])
    th_r = [sb.media_ponderada(d["pv1"], d[f"peso_r{r}"]) for r in range(1, 9)]
    esperado = np.sqrt(np.sum((np.array(th_r) - theta) ** 2) / (8 * 0.25))
    _, se = sb.media_pv_brr(d, "pv", **KW)
    assert se == pytest.approx(esperado)


def test_reg_recupera_coeficientes():
    rng = np.random.default_rng(1)
    n = 5000
    x = rng.normal(size=n)
    d = pd.DataFrame({"peso": rng.uniform(1, 10, n)})
    for i in range(1, 4):
        d[f"pv{i}"] = 400 + 25 * x + rng.normal(scale=20, size=n)
    for r in range(1, 9):
        d[f"peso_r{r}"] = d["peso"] * rng.choice([0.5, 1.5], size=n)
    X = pd.DataFrame({"const": 1.0, "x": x})
    out = sb.reg_pv_brr(d, X, "pv", **KW)
    assert out.loc["x", "coef"] == pytest.approx(25, abs=1.0)
    assert out.loc["const", "coef"] == pytest.approx(400, abs=1.0)
    assert (out["se"] > 0).all()


def test_reg_coef_igual_a_beta_pv():
    d = hacer_datos()
    X = pd.DataFrame({"const": 1.0, "x": d["x"]})
    out = sb.reg_pv_brr(d, X, "pv", **KW)
    solo = sb.beta_pv(X, d, "pv", n_pv=3)
    assert out["coef"].to_numpy() == pytest.approx(solo.to_numpy())


def test_reg_exige_alineacion():
    d = hacer_datos()
    X = pd.DataFrame({"const": 1.0, "x": d["x"]}).iloc[:-1]
    with pytest.raises(ValueError):
        sb.reg_pv_brr(d, X, "pv", **KW)


def test_r2_ponderado_ajuste_perfecto():
    x = np.linspace(0, 1, 50)
    X = np.column_stack([np.ones(50), x])
    assert sb.r2_ponderado(X, 3 + 2 * x, np.ones(50)) == pytest.approx(1.0)


def test_sd_pond_pesos_iguales():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    assert sb.sd_pond(x, np.ones(4)) == pytest.approx(np.std(x))


def test_dummies_con_ref():
    d = pd.DataFrame({"c": ["a", "b", "c", "b"]})
    out = sb.dummies_con_ref(d, "c", "b", "v")
    assert list(out.columns) == ["v_a", "v_c"]
    with pytest.raises(ValueError):
        sb.dummies_con_ref(d, "c", "zzz", "v")


def test_media_por_grupo_n():
    d = hacer_datos()
    d["g"] = np.where(np.arange(len(d)) % 2 == 0, "A", "B")
    out = sb.media_por_grupo(d, "pv", ["g"], **KW)
    assert out["n"].sum() == len(d)
    assert set(out["g"]) == {"A", "B"}


# ------------------------------------------------------------- datos reales
RUTA = os.environ.get("SABER_DATA")
sin_datos = pytest.mark.skipif(not RUTA or not Path(RUTA).exists(),
                               reason="Defina SABER_DATA con la ruta a la base")

# (grado, sexo): (media, error estándar) de lectura, guía ICFES sección 5.4
GUIA = {
    (3, "F"): (402.0325, 2.825653), (3, "M"): (398.2147, 3.097142),
    (5, "F"): (404.1364, 1.100846), (5, "M"): (396.3050, 1.630569),
    (9, "F"): (401.9243, 2.073854), (9, "M"): (398.0341, 2.450003),
}


@pytest.fixture(scope="module")
def base():
    p = Path(RUTA)
    if p.suffix == ".parquet":
        df = pd.read_parquet(p)
    else:
        df = pd.read_csv(p, sep=";", encoding="utf-8", low_memory=False)
    df.columns = df.columns.str.lower().str.strip()
    return df.dropna(subset=["sexo"])  # la guía excluye los 14 sin sexo


@sin_datos
@pytest.mark.parametrize("clave", list(GUIA))
def test_guia_icfes_lectura(base, clave):
    grado, sexo = clave
    sub = base[(base["grado"].astype(int) == grado) & (base["sexo"] == sexo)]
    est, se = sb.media_pv_brr(sub, "lectpv")
    media_g, se_g = GUIA[clave]
    assert est == pytest.approx(media_g, abs=5e-5)
    assert se == pytest.approx(se_g, abs=5e-5)


@sin_datos
def test_wls_real_contra_statsmodels(base):
    sm = pytest.importorskip("statsmodels.api")
    d = base[(base["grado"].astype(int) == 5)].dropna(subset=["ind_hse", "nse"])
    X = pd.get_dummies(d[["ind_hse", "nse", "sexo"]], columns=["nse", "sexo"],
                       drop_first=True, dtype=float)
    X.insert(0, "const", 1.0)
    mod = sm.WLS(d["lectpv1"].to_numpy(), X, weights=d["peso"].to_numpy()).fit()
    mio = sb.beta_wls(X, d["lectpv1"], d["peso"])
    assert mio == pytest.approx(mod.params.values, abs=1e-8)
