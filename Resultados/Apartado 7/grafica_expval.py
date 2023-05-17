import matplotlib.pyplot as plt

# Leer datos dl archivo
x = []
K = []
error_x = []
error_K = []
with open('expval_posicion_lambda0.9.txt', 'r') as archivo:
    for linea in archivo:
        datos = linea.split(' ')
        x.append(float(datos[0]))
        K.append(float(datos[1]))
        error_x.append(float(datos[2]))
        error_K.append(float(datos[3]))

# Representar los datos simulados
plt.errorbar(x, K, xerr=error_x, yerr=error_K, color='#0E4749')

# Nombres de los ejes
plt.xlabel('t')
plt.ylabel('Valor esperado de la posición')

# Configurar leyenda
plt.legend()

# Mostrar y guardar la gráfica
plt.show()
#plt.savefig('grafica_expval_posicion_0.9.png', dpi=300, bbox_inches='tight')