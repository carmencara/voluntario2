import matplotlib.pyplot as plt

# Leer datos del archivo
x = []
K = []
error_x = []
error_K = []
with open('expval_energia_lambda0.9.txt', 'r') as archivo:
    for linea in archivo:
        datos = linea.split(' ')
        x.append(float(datos[0]))
        K.append(float(datos[1]))
        error_x.append(float(datos[2]))
        error_K.append(float(datos[3]))

# Representar los datos simulados
plt.errorbar(x, K, xerr=error_x, yerr=error_K, color='#AA4465', ecolor='#F59FBA', linestyle=' ', marker='.')

# Nombres de los ejes
plt.xlabel('t')
plt.ylabel('Valor esperado de la energía total')
plt.xlim([0,2000])

# Mostrar y guardar la gráfica
#plt.show()
plt.savefig('grafica_expval_energia_0.9.png', dpi=300, bbox_inches='tight')