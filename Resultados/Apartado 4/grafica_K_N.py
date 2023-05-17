import matplotlib.pyplot as plt

# Leer datos de cada archivo
N1 = []
K1 = []
error_N1 = []
error_K1 = []
with open('datos_lambda0.3.txt', 'r') as archivo:
    # Saltar las dos líneas comentadas
    archivo.readline()
    archivo.readline()
    for linea in archivo:
        datos = linea.split('\t')
        N1.append(float(datos[0]))
        K1.append(float(datos[1]))
        error_N1.append(float(datos[2]))
        error_K1.append(float(datos[3]))
N2 = []
K2 = []
error_N2 = []
error_K2 = []
with open('datos_lambda0.9.txt', 'r') as archivo:
    # Saltar las dos líneas comentadas
    archivo.readline()
    archivo.readline()
    for linea in archivo:
        datos = linea.split('\t')
        N2.append(float(datos[0]))
        K2.append(float(datos[1]))
        error_N2.append(float(datos[2]))
        error_K2.append(float(datos[3]))
N3 = []
K3 = []
error_N3 = []
error_K3 = []
with open('datos_lambda4.txt', 'r') as archivo:
    # Saltar las dos líneas comentadas
    archivo.readline()
    archivo.readline()
    for linea in archivo:
        datos = linea.split('\t')
        N3.append(float(datos[0]))
        K3.append(float(datos[1]))
        error_N3.append(float(datos[2]))
        error_K3.append(float(datos[3]))

# Trazar las gráficas
plt.errorbar(N1, K1, xerr=error_N1, yerr=error_K1, color='#AA4465', label='$\lambda=0.3$')
plt.errorbar(N2, K2, xerr=error_N2, yerr=error_K2, color='#0E4749', label='$\lambda=0.9$')
plt.errorbar(N3, K3, xerr=error_N3, yerr=error_K3, color='#E55812', label='$\lambda=4$')

# Nombres de los ejes
plt.xlabel('N')
plt.ylabel('K')

# Configurar leyenda
plt.legend(loc=(0.75, 0.7))

# Mostrar y guardar la gráfica
#plt.show()
plt.savefig('grafica_K_N.png', dpi=300, bbox_inches='tight')
