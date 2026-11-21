# Hallazgos

Cada resultado se clasifica por nivel de certeza: **verificado** (sale directo de los datos o de una comparación), **inferido** (consistente con los datos, pero no demostrado) y **descartado** (hipótesis que se probó y no se sostuvo). Todo es asociación en datos transversales y autorreportados.

## 1. Validación de la implementación

| Referencia | Resultado |
|---|---|
| Guía ICFES, media de lectura por grado y sexo (6 grupos) | Idéntica en media y error estándar |
| `statsmodels` (`WLS`, un PV) | Diferencia máxima ≈ 2,7e-12 |
| Stata, modelo 4 (N = 67.321), 92 coeficientes (46 por área, constante incluida) | Diferencia máxima de 2e-05 en coeficientes y 0,0 (a 5 decimales) en errores estándar |
| `pytest` | 21 de 21 |

Las diferencias de 2e-05 en coeficientes son compatibles con que Stata importó los datos en precisión simple; es una explicación probable, no comprobada. La comparación completa queda en `reports/validacion_python_stata.csv`.

## 2. Calidad de datos (verificado)

- 108.074 registros; 14 sin `sexo` (de ahí los 108.060 de la tesis).
- Rango del puntaje 100-700 (la tesis dice 0-500). Índices con media ≈ 5 y DE ≈ 2; mínimo 0, valores 0 reales (rinden menos), no centinelas.
- 6 regiones (la tesis dice 5 en el texto).
- Nulos estructurales en 3º (100 %): `sub_se_8`, `ind_fe_4`, educación del cuidador (padre y madre).
- Pesos con cola larga (295 registros con peso > 200, 0,3 %). No se recortan.

## 3. Réplica de la Tabla 4

**Verificado**
- El modelo solo incluye 5º y 9º (N ≈ 67.300). La tesis describe tres grados y no reporta el N del modelo.
- La fila "Mujer" corresponde a hombres (+16 en matemáticas, -6 en lectura). La referencia de grado es 5º, no 3º.
- Los controles (sexo, sector, rural, `socioe`, `trabajo_equipo`, `entorno`, `pref_apr`, `invol_parental`, grado 9º, años previos) se reproducen dentro de 2 errores estándar.
- NSE se reproduce de forma aproximada (1,6 a 3,2 EE con `sub_fe_9`).

**Inferido**
- `doc_est` (11,26 / 13,50) y `percep_colegio` (-12,48 / -14,55) solo se reproducen añadiendo `sub_fe_9`: diferencia de 14 a 21 EE sin él, de 0,9 a 1,5 EE con él. `percep_colegio` es un compuesto que contiene a `doc_est` y a `sub_fe_9`; incluir compuesto y componente produce coeficientes de signos opuestos que no se pueden interpretar por separado. La tesis no declara esa variable.

| Variable | Tesis (mate / lect) | Réplica base | Con `sub_fe_9` |
|---|---|---|---|
| `doc_est` | 11,26 / 13,50 | 2,82 / 5,41 | 10,12 / 12,04 |
| `percep_colegio` | -12,48 / -14,55 | +0,65 / -1,96 | -10,81 / -12,38 |
| `socioe` | 4,32 / 3,76 | 4,28 / 3,71 | 4,22 / 3,68 |

**R²** (inferido): la tesis reporta 0,310 / 0,306. Con un PV por vez, 0,27 / 0,26. Con el promedio de los 10 PV como dependiente, 0,3076 / 0,2991. Es compatible con que la tesis usó el promedio, pero no lo demuestra. Ese procedimiento deja los coeficientes iguales y subestima los errores estándar (de 1 % en `publico` a 38 % en `sub_fe_9`; 28 % en `doc_est`).

**Descartado como explicación del R²:** faltaban carnes y años previos (solo lo subieron 0,01), la herramienta (Stata y Python coinciden), el cálculo del R².

## 4. Resultado propio: relación estudiante-docente

Modelo limpio (`doc_est` y `gusto_colegio` por separado, sin el compuesto):

| | Matemáticas | Lectura |
|---|---|---|
| Coeficiente (5º y 9º, N = 67.321) | 3,17 (EE 0,30) | 4,08 (EE 0,24) |
| Modelo ampliado, 3 grados (N = 91.133) | 2,96 (EE 0,30) | 3,94 (EE 0,20) |
| Efecto por 1 DE, por grado (3º / 5º / 9º) | 0,089 / 0,091 / 0,079 | 0,117 / 0,114 / 0,090 |

- Positivo y de tamaño pequeño (≈ 0,08 a 0,12 DE) en los tres grados. Estable por grado (con `socioe` libre por grado, las interacciones `doc_est × grado` tienen |t| < 2).
- Robusto a: excluir pesos extremos (cambios < 0,4 puntos), añadir carnes y años previos (< 0,15 puntos), y cambiar de especificación.
- El efecto es pequeño frente a la brecha socioeconómica (NS4 frente a NS1 ≈ 30 a 37 puntos).

## 5. Resultado propio: habilidades socioemocionales

**Verificado**
- `ind_hse` se compone de dimensiones distintas por grado: 3º, autoconciencia + autorregulación; 5º y 9º, automotivación + autorregulación. Solo autorregulación es común a los tres grados. La comparación directa del índice general entre grados no es válida.
- Coeficiente de `socioe` por grado (mate): 11,60 / 7,19 / 0,45. En 9º no hay asociación detectable.

**Con autorregulación (común a los tres grados)**, efecto por 1 DE: 0,134 / 0,053 / 0,020 (mate) y 0,144 / 0,073 / 0,025 (lectura). La asociación decrece, pero no se puede decir si es un cambio real con la edad o un cambio en cómo responden los estudiantes.

## 6. Hipótesis descartadas o retiradas

- **Heterogeneidad de `doc_est` por grado:** apareció con una especificación que forzaba `socioe` a ser igual en todos los grados; desapareció al liberarlo.
- **Los índices de clima aumentan/reducen los coeficientes de NSE:** NS3 y NS4 cambian ≤ 1 % al incluir o quitar los índices. El salto visto en una corrida de Stata venía de un error de código (se omitieron todas las variables numéricas).
- **Mitigación de brechas por el clima:** las seis interacciones `doc_est × NSE` tienen |t| < 2. Resultado nulo, con poca potencia en NS4 (7 % de la muestra) y sensible a la especificación. No respalda ni descarta la afirmación de la tesis.

## 7. Limitaciones

- Datos transversales; los dos índices los responde el mismo estudiante que presenta la prueba. Posible causalidad inversa (ej.: `invol_parental` es negativo).
- **No respuesta:** 6.890 registros sin ningún índice rinden unos 24 puntos menos (comparación simple, un PV, sin ponderar). La muestra del modelo supera a la población en 5,8 puntos (405,8 frente a 400,0). Los pesos no corrigen esto.
- Los IC de los efectos estandarizados tratan las DE como constantes y subestiman levemente la incertidumbre.
- No se evaluó el efecto de la estructura por colegio (efectos fijos o modelo multinivel), extensión natural.
- Sobre la tesis: lo verificado son problemas de reporte y de alcance. Que haya incluido `sub_fe_9` es una inferencia.
