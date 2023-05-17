import matplotlib.pyplot as plt
import numpy as np
import math
pi = 3.14159265359
ciclos = 50

# Leer datos dl archivo
x = []
K = []
error_x = []
error_K = []
with open('datos_N1000.txt', 'r') as archivo:
    # Saltar las dos líneas comentadas
    archivo.readline()
    archivo.readline()
    for linea in archivo:
        datos = linea.split('\t')
        x.append(float(datos[0]))
        K.append(float(datos[1]))
        error_x.append(float(datos[2]))
        error_K.append(float(datos[3]))

# Representar los datos simulados
plt.errorbar(x, K, xerr=error_x, yerr=error_K, color='#0E4749', linestyle='', marker='.', label='Datos simulados')

# Definir la curva teórica
def funcion(x):
    if x<1:
        return 4*(1-x)/(4.*(1-x)+(x*np.sin(2*pi*ciclos*np.sqrt(1-x)/5))**2)
    if x>1:
        return 4*(x-1)/(4.*(x-1)+(x*np.sinh(2*pi*ciclos*np.sqrt(x-1)/5))**2)
a = np.linspace(0,10,100)
y = [funcion(i) for i in a]

# Representar la curva teórica
plt.plot(a, y, label='Curva teórica', linestyle='--')

# Nombres de los ejes
plt.xlabel('$\lambda$')
plt.ylabel('K')

# Configurar leyenda
plt.legend()

# Mostrar y guardar la gráfica
plt.show()
plt.savefig('grafica_lambda_K.png', dpi=300, bbox_inches='tight')