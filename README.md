# Clima de aula y rendimiento en Saber 3°, 5° y 9° (2022): réplica y auditoría

Réplica en **Python y Stata** de la tesis *"Incidencia del Clima de Aula en el Rendimiento Escolar en Colombia"* (Maestría en Economía y Política de Educación, U. Externado, 2025), con inferencia correcta para encuestas complejas: **10 valores plausibles (PV) + 84 réplicas balanceadas (BRR, Fay = 0,5)**.

## Pregunta

¿Cómo se asocian la relación estudiante-docente y las habilidades socioemocionales con el puntaje en matemáticas y lectura, controlando por contexto? ¿Se sostienen los resultados de la tesis al replicarlos?

## Resultados en breve

1. **El modelo de la tesis solo usa 5º y 9º.** El grado 3 tiene 100 % de nulos en tres variables del modelo y desaparece por eliminación de casos. Además, la fila "Mujer" de la Tabla 4 corresponde a hombres.
2. **El coeficiente de la relación docente-estudiante (11,3 en matemáticas) solo se aproxima al incluir un índice compuesto junto con uno de sus componentes** (`sub_fe_9`), que la tesis no declara. Con los componentes por separado baja a ≈ 3 puntos (≈ 0,08 a 0,12 DE del puntaje). Esto es una inferencia, no una prueba.
3. **La relación estudiante-docente se asocia positivamente con el puntaje en los tres grados**, con efecto pequeño y estable. Las habilidades socioemocionales también, pero la asociación decrece con el grado y el índice no mide lo mismo en todos los grados.
4. **El R² de la tesis (0,31) se reproduce regresando sobre el promedio de los 10 PV**, procedimiento que subestima los errores estándar entre 1 % y 38 %.
5. **No hay evidencia de que el clima de aula mitigue las brechas por nivel socioeconómico** (resultado nulo, con poca potencia en NS4).

Todo es **asociación**, no causalidad. Detalle, niveles de certeza y límites en [`reports/hallazgos.md`](reports/hallazgos.md).

## Validación de la implementación

| Referencia | Resultado |
|---|---|
| Guía metodológica ICFES (media de lectura por grado y sexo, 6 grupos) | Media y error estándar idénticos (tolerancia 5e-5) |
| `statsmodels` (`WLS`) | Coeficientes idénticos (diferencia ≈ 1e-12) |
| Stata (`svy` + `brrweight`, `mse`) | Mismo N (67.321); coeficientes y errores estándar del modelo 4 idénticos (dif. máx. 2e-05) |
| Tests (`pytest`) | 21 de 21 pasan |

## Estructura

```
saber-clima-aula/
├── README.md
├── requirements.txt
├── src/saber_brr.py          # media_pv_brr, reg_pv_brr, beta_wls, dummies_con_ref, ...
├── tests/test_saber_brr.py   # 14 sintéticos + 7 con datos reales
├── notebooks/analisis_replica.ipynb   # auditoría, réplica y extensión, de punta a punta
├── stata/replica.do          # modelo 4 completo, con asserts
├── reports/hallazgos.md
└── data/                     # vacío: ver "Datos"
```

## Datos

La base **no se incluye**. Descárguela de [DataIcfes](https://www.icfes.gov.co/investigaciones/data-icfes/) (Pruebas Saber 3°, 5°, 7°, 9°, aplicación 2022) y consulte su diccionario y guía metodológica. Fuente a citar: Icfes, DataIcfes, Repositorio de Datos Abiertos, Pruebas Saber 3°, 5°, 7°, 9°. El ICFES no se responsabiliza por los resultados obtenidos a partir de sus bases.

## Cómo reproducir

**Python**
```bash
pip install -r requirements.txt
export SABER_DATA=/ruta/a/Prueba_Saber_3579_2022.txt   # o un .parquet
pytest tests -v        # sin SABER_DATA, los 7 tests con datos reales se omiten
```
```python
import sys; sys.path.insert(0, 'src')
import saber_brr as sb
est, se = sb.media_pv_brr(df_subgrupo, 'lectpv')          # media con PV y BRR
res = sb.reg_pv_brr(df, X, 'matepv')                       # coef, se, t
```

**Stata**: edite las dos rutas `global` al inicio de `stata/replica.do` y ejecútelo **completo** (Ctrl+D). Genera un log y `coef_modelo4.csv`.

## Método (guía ICFES, ecuaciones 1 a 5)

- Estimación: promedio de las 10 estimaciones, una por PV.
- Varianza de muestreo: BRR con corrección de Fay, centrada en la estimación principal (equivale a `mse` en Stata), promediada sobre PV.
- Varianza de imputación: varianza entre las 10 estimaciones.
- Total: `V_muestreo + (1 + 1/M) · V_imputación`.
- Las regresiones WLS se resuelven con álgebra matricial (850 regresiones por modelo en segundos).

## Uso de inteligencia artificial

Este proyecto se desarrolló con la asistencia de **Claude** (Anthropic), usado como asistente conversacional durante noviembre de 2026.

**Qué hizo la IA**
- Revisó los documentos de contexto (tesis, guía metodológica y diccionario del ICFES) y propuso el plan de trabajo por fases.
- Redactó borradores de código en Python y Stata, y de los tests, el do-file, el cuaderno y la documentación.
- Ayudó a interpretar resultados, a plantear hipótesis y a señalar cautelas (por ejemplo, la redundancia entre índices compuestos y sus componentes).

**Qué hizo el autor**
- Ejecutó todo el código sobre los datos reales, en Colab y en Stata, y revisó cada salida.
- Decidió el alcance, las especificaciones y qué conclusiones entran al informe.
- Es responsable del contenido final.

**Cómo se controló el riesgo de errores**
- Ningún resultado se aceptó sin contrastarlo con una fuente externa: las medias y errores estándar de la guía del ICFES, `statsmodels` y una implementación independiente en Stata (ver "Validación de la implementación").
- Las afirmaciones se clasifican por certeza (verificado, inferido, descartado) en [`reports/hallazgos.md`](reports/hallazgos.md).
- Varios errores del asistente se detectaron al ejecutar y se corrigieron: instrucciones de Stata que dejaban variables vacías al ejecutar por bloques, una cifra de conteo mal calculada y conclusiones sobre heterogeneidad por grado que luego se retiraron al contrastarlas con los datos.

**Qué no hizo la IA:** no tuvo acceso a la base de datos. Trabajó sobre capturas de pantalla y texto que el autor le compartió, por lo que no pudo ejecutar el código con los datos reales.

## Advertencias

- Los pesos no se recortan: el ICFES calibró el peso y las réplicas conjuntamente.
- Los coeficientes de la regresión sobre el promedio de PV coinciden con los de la combinación, pero sus errores estándar no son válidos.
