Este repositorio contiene los archivos relacionados con el problema voluntario de Física Computacional, correspondiente al tema 3
(resolución de la ecuación de Schrödinger). Los objetivos de este programa son:

1. Estima el valor del coeficiente de transmisión K = mT/m, d ́onde mT es el número de veces que hemos detectado
la partícula a la derecha del potencial. Para ello simularemos el sistema al menos 10E3 veces.
2. Compara el valor de K = mT/m obtenido y el de la probabilidad pD(nD) en el punto máximo obtenido t=nD (promediada sobre el número
de simulaciones). Explica la relación observada entre estos números con argumentos físicos sencillos (clásicos o cuánticos). Explica
por qué nuestro método de maximización tiene sentido.
3. Estima el error estadístico de nuestro estimador K. Pista: puedes usar el teorema de Hoeffding para variables de Bernouilli, o los
métodos de cálculo de errores que vimos en el tema de algoritmos Monte Carlo.
4. Estudia la dependencia en N de K realizando simulaciones para N = 500, 1000, 2000.
5. Estudiar la dependencia en V (x) de K, usando valores λ = 0,1; 0,3; 0,5; 1; 5; 10.
6. Comparar con resultados teóricos (puedes usar simplificaciones como soluciones analíticas para funciones de onda más sencillas y
comparar los órdenes de magnitud).
7. Calcula los valores esperados de distintos observables (posición, momento, energía cinética, energía total, ...) con errores y
discutir los resultados.
8. Estudia el comportamiento de K en función del número de barreras n y cuantifica dicho comportamiento realizando  un ajuste con
una función de prueba adecuada. Pista: puedes pintar primero la función gráficamente y usa dicha información para determinar si puede
ser apropiado realizar un ajuste polinomial, logarítmico o exponencial; tambi ́en se pueden usar métricas para comprobar la bondad
del ajuste.

Se han realizado dos programas en C++:
- transmision.cpp (para los apartados 1-7)
- barreras.cpp (para el apartado 8)
En ambos se ha necesitado la librería GSL para la generación de números aleatorios (gsl_rng.h).

Los ficheros con los datos que se extraen de las simulaciones se encuentran en la carpeta "Resultados", que está organizada en las
distintas tareas que se requieren en este problema voluntario:
- Apartados 1, 2 y 3: "potencial.txt" es la forma del potencial, que se grafica con "grafica_potencial.py". Además, los datos relativos
  al coeficiente de transmisión están en "coeficiente_transmision.txt".
- Apartado 4: los archivos "coeficiente_transmision_lambda%N*.txt" salen de las simulaciones directamente. En "datos_lambda%.txt" se
  resumen los datos necesarios para graficar K en función de N.
- Apartados 5 y 6: los archivos "coeficiente_transmision_lambda%.txt" salen de las simulaciones directamente. En "datos_N1000.txt"
  están los valores de K en función de lambda.
- Apartado 7: están los resultados para los valores esperados que salen de las simulaciones. Con "grafica_expval.py" se grafican estos
  valores junto a su error.
- Apartado 8: los archivos "X_coeficiente_transmision.txt" son los resultados relativos al coeficiente de transmisión para X barreras
  (X=1,2...10). Para representar K en función del número de barreras se usa "grafica_barreras_K.py" y el ajuste se hace con Gnuplot,
  en concreto con "grafica_barreras_K.plt" ("fit.log" son los datos relativos al ajuste). Además, para representar la forma del
  potencial "barreras_potencial.txt" se necesita "grafica_potencial.py".
  
