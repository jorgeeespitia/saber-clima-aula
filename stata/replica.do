*==============================================================================
* replica.do  -  Replicación en Stata (Saber 3°, 5°, 7° y 9°, 2022)
* Modelo 4: sub_fe_9 + años previos + carnes, 5º y 9º, 10 PV + BRR (Fay 0.5, mse)
* Ejecutar COMPLETO (Ctrl+D). Los local no persisten entre ejecuciones parciales.
*==============================================================================
clear all
set more off

*--- 0. Rutas (editar solo aquí) ----------------------------------------------
global datos "\ruta\a\datos"
global salida "C:\ruta\a\guardar"
capture log close
log using "$salida\replica_log.txt", text replace

*--- 1. Carga y limpieza -------------------------------------------------------
import delimited "$datos", delimiter(";") encoding(UTF-8) varnames(1) clear

* "NA" llega como texto: se convierte en faltante
foreach v of varlist _all {
    capture confirm string variable `v'
    if _rc==0 {
        quietly replace `v' = "" if `v' == "NA"
    }
}

* Índices llegan como texto: se convierten a numérico (sin force, para no perder datos)
destring ind_hse sub_fe_11 sub_fe_9 ind_fe_3 ind_fe_4 ind_fe_5 ind_fe_6 sub_se_8 ind_se_3, replace

assert _N == 108074
count if missing(ind_hse)
assert r(N) == 7727

*--- 2. Diseño muestral --------------------------------------------------------
svyset [pw=peso], brrweight(peso_r1-peso_r84) fay(.5) vce(brr) mse

*--- 3. Validación contra la guía ICFES (media de lectura, 6 grupos) -----------
local gs   3 3 5 5 9 9
local ss   F M F M F M
local mus  402.0325 398.2147 404.1364 396.3050 401.9243 398.0341
local sds  2.825653 3.097142 1.100846 1.630569 2.073854 2.450003

forvalues i = 1/6 {
    local g  : word `i' of `gs'
    local s  : word `i' of `ss'
    local mu : word `i' of `mus'
    local sd : word `i' of `sds'
    forvalues m = 1/10 {
        quietly svy, subpop(if grado==`g' & sexo=="`s'"): mean lectpv`m'
        matrix b`m' = e(b)
        matrix V`m' = e(V)
    }
    scalar bbar = 0
    scalar Vbar = 0
    forvalues m = 1/10 {
        scalar bbar = bbar + b`m'[1,1]/10
        scalar Vbar = Vbar + V`m'[1,1]/10
    }
    scalar B = 0
    forvalues m = 1/10 {
        scalar B = B + (b`m'[1,1] - bbar)^2/9
    }
    scalar se = sqrt(Vbar + (1 + 1/10)*B)
    display "Grado `g' `s': media " %9.4f bbar "  EE " %9.4f se "  | guía: `mu'  `sd'"
    assert abs(bbar - `mu') < 1e-4 & abs(se - `sd') < 1e-4
}
display "VALIDACIÓN GUÍA ICFES: OK"

*--- 4. Variables del modelo ---------------------------------------------------
gen hombre        = (sexo=="M") if sexo!=""
gen publico       = (naturaleza_establecimiento=="PUBL")
gen publico_rural = (seczon=="Público - Rural")
gen grado9        = (grado==9)
gen internet_si   = (i1775476=="Sí") if i1775476!=""

* Categoría "NA" explícita solo en años previos y carnes
replace i1817204 = "NA" if missing(i1817204)
replace i1775443 = "NA" if missing(i1775443)

* encode ordena alfabéticamente; la categoría 1 es la referencia en todas
* (NS1, Región Caribe, Entre 1 y 2 horas, 0 a 10 libros, Bachillerato completo, 0 años, 0 días)
foreach v in nse region i1775328 i1775434 i1864393 i1817204 i1775443 {
    encode `v', gen(`v'_n)
}

*--- 5. Muestra ----------------------------------------------------------------
local xnum sub_fe_11 ind_hse grado9 hombre publico publico_rural ///
    ind_fe_3 ind_fe_4 ind_fe_5 ind_fe_6 sub_se_8 internet_si sub_fe_9
local xfac i.nse_n i.region_n i.i1775328_n i.i1775434_n i.i1864393_n i.i1817204_n i.i1775443_n

egen nmiss = rowmiss(`xnum' nse_n region_n i1775328_n i1775434_n i1864393_n)
gen byte touse = inlist(grado,5,9) & nmiss==0
count if touse
assert r(N) == 67321
display "xnum = `xnum'"

*--- 6. Regresión con 10 PV y BRR (ecuaciones 1-5 de la guía) -------------------
tempname h
postfile `h' str5 area str40 termino coef se using "$salida\coef_modelo4.dta", replace

foreach area in mate lect {
    forvalues m = 1/10 {
        quietly svy, subpop(touse): regress `area'pv`m' `xnum' `xfac'
        matrix b`m' = e(b)
        matrix V`m' = e(V)
        scalar r2_`area'`m' = e(r2)
    }
    matrix bbar = b1*0
    matrix Vbar = V1*0
    forvalues m = 1/10 {
        matrix bbar = bbar + b`m'/10
        matrix Vbar = Vbar + V`m'/10
    }
    matrix Bm = V1*0
    forvalues m = 1/10 {
        matrix dd = b`m' - bbar
        matrix Bm = Bm + dd'*dd/9
    }
    matrix T = Vbar + (1 + 1/10)*Bm

    local nombres : colnames bbar
    local k = colsof(bbar)
    display _n "Resultados `area'"
    display "{hline 45}"
    forvalues j = 1/`k' {
        local nom : word `j' of `nombres'
        display %-22s "`nom'" %9.3f bbar[1,`j'] %9.3f sqrt(T[`j',`j'])
        post `h' ("`area'") ("`nom'") (bbar[1,`j']) (sqrt(T[`j',`j']))
    }

    scalar r2m = 0
    forvalues m = 1/10 {
        scalar r2m = r2m + r2_`area'`m'/10
    }
    display "R2 promedio de PV individuales (`area'): " %6.3f r2m
}
postclose `h'

*--- 7. R2 con el promedio de los PV (reproduce el R2 de la tesis) --------------
capture drop mate_prom lect_prom
egen mate_prom = rowmean(matepv1-matepv10)
egen lect_prom = rowmean(lectpv1-lectpv10)

svy, subpop(touse): regress mate_prom `xnum' `xfac'
display "R2 mate (promedio PV): " %6.4f e(r2)
svy, subpop(touse): regress lect_prom `xnum' `xfac'
display "R2 lect (promedio PV): " %6.4f e(r2)
* Nota: los coeficientes coinciden con la combinación de PV, pero los errores
* estándar son menores (subestimados 1 %-38 %) porque no incluyen la varianza de imputación.

*--- 8. Exportar coeficientes a CSV ----------------------------------------------
use "$salida\coef_modelo4.dta", clear
export delimited using "$salida\coef_modelo4.csv", replace

log close
