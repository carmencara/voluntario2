import matplotlib.pyplot as plt
import numpy as np

# Función exponencial
def f(x, a, b):
    return a * np.exp(-b * x)

# Leer datos del archivo
x = []
K = []
error_x = []
error_K = []
with open('barreras_K.txt', 'r') as archivo:
    archivo.readline()  # Saltar una línea comentada
    for linea in archivo:
        datos = linea.split('\t')
        x.append(float(datos[0]))
        K.append(float(datos[1]))
        error_x.append(float(datos[2]))
        error_K.append(float(datos[3]))

# Definir los valores de a y b
a = 0.764418
b = 0.192254

# Generar puntos para la curva ajustada
x_fit = np.linspace(min(x), max(x), 100)
y_fit = f(x_fit, a, b)

# Representar los datos simulados
plt.errorbar(x, K, xerr=error_x, yerr=error_K, color='#0E4749', linestyle='', marker='.', label='Datos simulados')

# Representar la curva ajustada
plt.plot(x_fit, y_fit, color='#14BDEB', label='Ajuste exponencial')

# Nombres de los ejes
plt.xlabel('n')
plt.ylabel('K')

# Configurar leyenda
plt.legend()

# Mostrar y guardar la gráfica
#plt.show()
plt.savefig('grafica_barreras_K.png', dpi=300, bbox_inches='tight')
