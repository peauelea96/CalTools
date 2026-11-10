# Validación técnica de CalTools 0.4.0

## Cambios

- Preparación dispersa común para controles, intervalos y cotas individuales.
- OSQP con Hessiana diagonal, estados precisos y warm start primal.
- Dykstra ponderado con correcciones compactas y comprobación KKT.
- Newton dual semisuave experimental, regularización, Armijo y respaldo OSQP opcional.
- Selección automática documentada; rutas históricas conservadas para entradas pequeñas.
- calibration_diagnostics: totales, violaciones, lambda, pesos de entrada/salida y g.
- DESCRIPTION, NAMESPACE, ayuda roxygen, ejemplos y tests actualizados.

## Validación

Los tests históricos y nuevos suman 793 comprobaciones (478 históricas y 315
nuevas). Sin fallos ni advertencias. Los nuevos casos cubren comparaciones con
quadprog, igualdades, desigualdades activas/inactivas, intervalos, cotas fijas y
vectoriales, redundancia, ceros, totales negativos, dispersión, no convergencia,
infactibilidad, q no constante y un caso mal condicionado. En el caso mal
condicionado se verifica la honestidad del estado; no se exige convergencia.

R CMD check --no-manual --no-build-vignettes: sin errores ni advertencias.
Una NOTE: rmarkdown y sampling no están disponibles en el entorno. No se
compiló el manual PDF en esta comprobación. R 4.3.3, Matrix 1.6.5, osqp 1.0.0.
Los registros se incluyen en validation/.

### Resultados medidos en este entorno

Medianas de tiempo total (segundos), tres repeticiones, perfil de 16 controles,
6 desigualdades, densidad .1 y cotas estrechas. Tolerancia 1e-6 y máximo 1.000
iteraciones por método. El asterisco indica una solución aproximada de OSQP,
por lo que su tiempo no representa haber alcanzado la misma precisión.

| n | OSQP | Dykstra | Newton dual |
|---:|---:|---:|---:|
| 1000 | 0.007 | 0.014 | 0.007 |
| 3000 | 0.015 | 0.015 | 0.010 |
| 10000 | 0.114 | 0.058 | 0.022 |
| 50000 | 1.806* | 0.236 | 0.098 |
| 100000 | 4.395* | 0.444 | 0.226 |

De 144 ejecuciones realizadas, 138 pasaron KKT y 6 quedaron como aproximadas;
otras 36 combinaciones de quadprog se omitieron por el límite de tamaño.
En el perfil estrecho, el máximo residuo KKT de OSQP fue 0.00539 para
50.000 unidades y 0.00964 para 100.000. Los CSV conservan objetivos, memoria
asignada por R, cotas activas y tiempos separados de todas las ejecuciones.

Estas muestras favorecen estructuras con pocos controles y no demuestran
superioridad universal. La validación de Newton sigue siendo limitada; faltan
ensayos extensos con datos reales, conjuntos casi infactibles y escalas extremas.

## Límites pendientes

Newton es experimental; no hay garantía global. El control-space Gram denso
limita Dykstra/Newton cuando hay muchos controles. La pseudoinversa de Dykstra
usa tolerancia de rango fija 1e-12; para degeneración extrema se debe revisar
el estado y reescalar controles. OSQP 0.6 tiene adaptador de compatibilidad,
pero las pruebas de esta entrega se ejecutaron con 1.0.0. Warm starts entre
llamadas no reutilizan factorizaciones. Los certificados de infactibilidad son
numéricos, no exactos. Dykstra y Newton no prueban infactibilidad.

La memoria medida es asignación acumulada de R, no pico RSS ni memoria nativa
completa de OSQP. Los benchmarks tienen solo tres perfiles sintéticos y tres
repeticiones en un único entorno. Quadprog se limita a n<=2000, de modo que
no se ha forzado una comparación densa a 100.000 unidades.

La ruta legacy mantiene errores ante no convergencia y su status original;
las rutas nuevas devuelven estados explícitos. Las notas de compatibilidad
del README identifican este cambio para auto en matrices dispersas o n>2000.
