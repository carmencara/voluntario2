import matplotlib.pyplot as plt

N = 1000
nbarreras = 3
ancho=N+2*N/5*(nbarreras-1)

# Leer datos del archivo
x_valores = []
v_valores = []
with open('barreras_potencial.txt', 'r') as archivo:
    for linea in archivo:
        x, v = linea.strip().split()  # Suponiendo que las columnas están separadas por espacios
        x_valores.append(float(x))
        v_valores.append(float(v))

# Crear la figura y los ejes
fig, ax = plt.subplots()

# Trazar la gráfica
ax.plot(x_valores, v_valores, color='#AA4465')

# Trazar las líneas verticales
ax.axvline(x=0, color='#AA4465')
ax.axvline(x=ancho, color='#AA4465')

# Nombres de los ejes
ax.set_xlabel('x')
ax.set_ylabel('V(x)')
ax.set_ylim([-0.01,1.3])

# Mostrar y guardar la gráfica
#plt.show()
plt.savefig('grafica_barreras_potencial.png', dpi=300, bbox_inches='tight')
