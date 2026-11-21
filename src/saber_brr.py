"""
saber_brr: inferencia para encuestas complejas con valores plausibles (PV)
y réplicas balanceadas (BRR) para las Pruebas Saber 3°, 5°, 7° y 9° (ICFES).

Sigue la guía metodológica del ICFES (ecuaciones 1 a 5):
  - Estimación puntual: promedio de las M estimaciones (una por PV).
  - Varianza de muestreo: BRR con corrección de Fay, centrada en la
    estimación principal (equivale a ``mse`` en Stata), promediada sobre PV.
  - Varianza de imputación: varianza entre las M estimaciones.
  - Varianza total: V_muestreo + (1 + 1/M) * V_imputación.

Valores por defecto: M = 10 PV, G = 84 réplicas, Fay = 0.5.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

M_PV = 10
G_REP = 84
FAY = 0.5


# ---------------------------------------------------------------- utilidades
def _cols_pv(prefix: str, n_pv: int) -> list[str]:
    return [f"{prefix}{i}" for i in range(1, n_pv + 1)]


def _cols_rep(prefix: str, n_rep: int) -> list[str]:
    return [f"{prefix}{i}" for i in range(1, n_rep + 1)]


def _combinar_varianza(var_brr, estimaciones, n_pv, axis=0):
    """Varianza total = promedio(BRR) + (1 + 1/M) * varianza entre PV."""
    return var_brr.mean(axis=axis) + (1 + 1 / n_pv) * estimaciones.var(axis=axis, ddof=1)


# -------------------------------------------------------------------- medias
def media_ponderada(y, w) -> float:
    """Media ponderada de ``y`` con pesos ``w``."""
    y, w = np.asarray(y, dtype=float), np.asarray(w, dtype=float)
    return float(np.sum(w * y) / np.sum(w))


def media_pv_brr(d: pd.DataFrame, pv_prefix: str, peso: str = "peso",
                 rep_prefix: str = "peso_r", n_pv: int = M_PV,
                 n_rep: int = G_REP, fay: float = FAY) -> tuple[float, float]:
    """Media con valores plausibles y error estándar total (BRR + imputación).

    Parameters
    ----------
    d : DataFrame con las columnas de PV, peso y réplicas.
    pv_prefix : prefijo de los PV, p. ej. ``'lectpv'`` (lectpv1 ... lectpv10).

    Returns
    -------
    (estimación, error estándar)
    """
    w = d[peso].to_numpy(dtype=float)
    W_rep = d[_cols_rep(rep_prefix, n_rep)].to_numpy(dtype=float)  # n x G

    medias = np.empty(n_pv)
    var_brr = np.empty(n_pv)
    for m, pv in enumerate(_cols_pv(pv_prefix, n_pv)):
        y = d[pv].to_numpy(dtype=float)
        theta = media_ponderada(y, w)
        medias[m] = theta
        theta_r = (W_rep * y[:, None]).sum(axis=0) / W_rep.sum(axis=0)
        var_brr[m] = np.sum((theta_r - theta) ** 2) / (n_rep * (1 - fay) ** 2)

    var_total = _combinar_varianza(var_brr, medias, n_pv)
    return float(medias.mean()), float(np.sqrt(var_total))


def media_por_grupo(d: pd.DataFrame, pv_prefix: str, by: list[str],
                    **kwargs) -> pd.DataFrame:
    """Aplica ``media_pv_brr`` a cada grupo definido en ``by``."""
    filas = []
    for clave, sub in d.groupby(by, observed=True):
        clave = clave if isinstance(clave, tuple) else (clave,)
        est, se = media_pv_brr(sub, pv_prefix, **kwargs)
        filas.append((*clave, est, se, len(sub)))
    return pd.DataFrame(filas, columns=[*by, "media", "se", "n"])


# ----------------------------------------------------------------- regresión
def beta_wls(X, y, w) -> np.ndarray:
    """Coeficientes WLS por solución cerrada: (X'WX)^-1 X'Wy."""
    X = np.asarray(X, dtype=float)
    XtW = X.T * np.asarray(w, dtype=float)  # k x n
    return np.linalg.solve(XtW @ X, XtW @ np.asarray(y, dtype=float))


def reg_pv_brr(d: pd.DataFrame, X: pd.DataFrame, pv_prefix: str,
               peso: str = "peso", rep_prefix: str = "peso_r",
               n_pv: int = M_PV, n_rep: int = G_REP,
               fay: float = FAY) -> pd.DataFrame:
    """Regresión WLS con PV y BRR. Devuelve coef, se y t por columna de ``X``.

    ``X`` debe estar alineada con ``d`` (mismas filas, mismo orden) e incluir
    la constante. Hace n_pv * (1 + n_rep) regresiones (850 por defecto).
    """
    if len(X) != len(d):
        raise ValueError("X y d deben tener el mismo número de filas")
    Xm = X.to_numpy(dtype=float)
    w = d[peso].to_numpy(dtype=float)
    W_rep = d[_cols_rep(rep_prefix, n_rep)].to_numpy(dtype=float)
    k = Xm.shape[1]

    betas = np.empty((n_pv, k))
    var_brr = np.empty((n_pv, k))
    for m, pv in enumerate(_cols_pv(pv_prefix, n_pv)):
        y = d[pv].to_numpy(dtype=float)
        b = beta_wls(Xm, y, w)
        betas[m] = b
        b_rep = np.array([beta_wls(Xm, y, W_rep[:, r]) for r in range(n_rep)])
        var_brr[m] = ((b_rep - b) ** 2).sum(axis=0) / (n_rep * (1 - fay) ** 2)

    est = betas.mean(axis=0)
    se = np.sqrt(_combinar_varianza(var_brr, betas, n_pv))
    return pd.DataFrame({"coef": est, "se": se, "t": est / se}, index=X.columns)


def beta_pv(X: pd.DataFrame, d: pd.DataFrame, pv_prefix: str,
            peso: str = "peso", n_pv: int = M_PV) -> pd.Series:
    """Solo coeficientes (promedio de los PV), sin errores estándar. Es rápido."""
    Xm, w = X.to_numpy(dtype=float), d[peso].to_numpy(dtype=float)
    b = np.mean([beta_wls(Xm, d[c].to_numpy(dtype=float), w)
                 for c in _cols_pv(pv_prefix, n_pv)], axis=0)
    return pd.Series(b, index=X.columns)


def r2_ponderado(X, y, w) -> float:
    """R² ponderado de una regresión WLS."""
    X, y, w = (np.asarray(a, dtype=float) for a in (X, y, w))
    res = y - X @ beta_wls(X, y, w)
    ybar = np.sum(w * y) / np.sum(w)
    return float(1 - np.sum(w * res ** 2) / np.sum(w * (y - ybar) ** 2))


# --------------------------------------------------------- desviaciones, dummies
def sd_pond(x, w) -> float:
    """Desviación estándar ponderada (poblacional)."""
    x, w = np.asarray(x, dtype=float), np.asarray(w, dtype=float)
    m = np.sum(w * x) / np.sum(w)
    return float(np.sqrt(np.sum(w * (x - m) ** 2) / np.sum(w)))


def sd_pv(d: pd.DataFrame, pv_prefix: str, peso: str = "peso",
          n_pv: int = M_PV) -> float:
    """DE ponderada del puntaje, promediada sobre los PV."""
    w = d[peso].to_numpy(dtype=float)
    return float(np.mean([sd_pond(d[c], w) for c in _cols_pv(pv_prefix, n_pv)]))


def dummies_con_ref(d: pd.DataFrame, col: str, ref: str, prefijo: str) -> pd.DataFrame:
    """Dummies de ``col`` eliminando explícitamente la categoría de referencia.

    Se evita ``drop_first=True`` porque elegiría la primera categoría
    alfabética, que puede no ser la referencia deseada.
    """
    dum = pd.get_dummies(d[col].astype(str), prefix=prefijo, dtype=float)
    ref_col = f"{prefijo}_{ref}"
    if ref_col not in dum.columns:
        raise ValueError(f"No existe la referencia {ref_col!r}; "
                         f"categorías: {list(dum.columns)}")
    return dum.drop(columns=ref_col)
