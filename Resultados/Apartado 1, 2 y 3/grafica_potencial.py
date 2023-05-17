import matplotlib.pyplot as plt

# Leer datos del archivo
x_valores = []
v_valores = []
with open('potencial.txt', 'r') as archivo:
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
ax.axvline(x=1000, color='#AA4465')

# Nombres de los ejes
ax.set_xlabel('x')
ax.set_ylabel('V(x)')
ax.set_ylim([-0.001,0.09])

# Mostrar y guardar la gráfica
#plt.show()
plt.savefig('grafica_potencial.png', dpi=300, bbox_inches='tight')
