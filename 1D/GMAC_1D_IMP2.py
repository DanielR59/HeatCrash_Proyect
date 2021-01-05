import numpy as np
import matplotlib.pyplot as plt
import time
from matplotlib.animation import FuncAnimation
from Funciones1D import hdf5
import matplotlib.image as mpimg


def Laplaciano_Dirichlet(N,r,k):
    A = np.zeros((N,N))
    A[0,0] = 1 + r * (k[1]+k[2]) 
    A[0,1] = -r * k[2]
    for i in range(1,N-1):
        A[i,i] = 1 + r * (k[i+1]+k[i+2])
        A[i,i+1] = -r * k[i+2]
        A[i,i-1] = -r * k[i+1]
    A[N-1,N-2] = -r * k[N]; 
    A[N-1,N-1] = 1 + r * (k[N]+k[N+1]) 
    return A

def Laplaciano_Neumman(N,r,k):
    A = np.zeros((N,N))
    A[0,0] = 1 + r * (k[1]+k[2]) 
    A[0,1] = -r * k[2]
    for i in range(1,N-1):
        A[i,i] = 1 + r * (k[i+1]+k[i+2])
        A[i,i+1] = -r * k[i+2]
        A[i,i-1] = -r * k[i+1]
    A[N-1,N-2] = -r * k[N]; 
    A[N-1,N-1] = 1 + r * k[N]
    return A

#Definición de tipo de fronteras:
def tipo_Fronteras():
    tipoFS = int(input('''
        Elije el tipo de fronteras del sistema
        1 - Dirichlet
        2 - Dirichlet + Neumman
        >> '''))
    if tipoFS == 1:
        tipoFrontera = 'Dirichlet'
    elif tipoFS == 2:
        tipoFrontera = 'Neumman'
    else:
        tipoFrontera = 'Dirichlet'
    return tipoFrontera

#Definicion de como leer las variables:
def lectura_Variables(tipoFrontera):
	opcionVariables = int(input('''
		Escriba el numero de la opcion que prefiera:

	    1 - Leer las variables desde un archivo HDF5
	    2 - Ingresar las variables desde la terminal
	    3 - Usar variables predefinidas
	    >> '''))

	if opcionVariables == 1:
	    paramsArchivo = hdf5.leerParametros('Datos','tolerancia','a','b','k', 'Tmax', 'N','boundA', 'boundB')
	    tolerancia = paramsArchivo['tolerancia']
	    a = paramsArchivo['a']
	    b = paramsArchivo['b']
	    Tmax = paramsArchivo['Tmax']
	    N = paramsArchivo['N']
	    boundA = paramsArchivo['boundA']
	    boundB = paramsArchivo['boundB']
	#Caso variables desde la terminal:
	elif opcionVariables == 2:
	    tolerancia = float(input('Tolerancia del error > '))
	    a = float(input('Inicio de la barra > '))
	    b = float(input('Final de la barra > '))
	    N = int(input('Número de nodos > '))
	    Tmax = float(input('Tiempo total de simulación > '))
	    boundA = float(input('Condición Dirichlet en (a) > '))
	    boundB = float(input('Condición '+ tipoFrontera+' en (b)> '))
	#Caso variables predefinidas:
	else:
	    tolerancia = 1e-6 #Tolerancia del error
	    a = 0.0 #Inicio de la barra
	    b = 1.0 #fin de la barra
	    N = 50 #Numero de nodos
	    Tmax = 1.0 #Tiempo total
	    boundA = -1 #Condicion de frontera en el extremo 1
	    boundB = 1 #Condicion de frontera en el extremo 2

	return tolerancia, a, b, N, Tmax, boundA, boundB


def lectura_K(N):
	tipoK = int(input('''
		La difusividad térmica es variable?
		1 - Si
		2 - No
		>> '''))
	if tipoK == 1:
		tipoVariabilidad = int(input('''
			Selecciona la forma en la que varía k:
			1 - En un solo rango de puntos
			2 - En varios rangos de puntos
			3 - De forma aleatoria
			>> '''))
		if tipoVariabilidad == 1:
			k = np.ones(N+2)
			knodo1 = int(input('''
				(Se preseleccionó un valor de k contante igual a 1.)
				¿Desde qué nodo quieres cambiar el valor de k?
				>> '''))
			knodo2 = int(input('''
				...Hasta qué nodo quieres cambiar el valor de k?
				>> '''))
			valorNuevoK = int(input('''
				...¿Qué valor quieres definir en los nodos,'''+str(knodo1)+'''a'''+str(knodo2)+'''?
				>> ''')) 
			k[knodo1:knodo2] = valorNuevoK
		elif tipoVariabilidad == 2:
			k = np.ones(N+2)
			while True:
				knodo1 = int(input('''
				(Se preseleccionó un valor de k contante igual a 1.)
				¿Desde qué nodo quieres cambiar el valor de k?
				>> '''))
				knodo2 = int(input('''
				...Hasta qué nodo quieres cambiar el valor de k?
				>> '''))
				valorNuevoK = int(input('''
				...¿Qué valor quieres definir en los nodos,'''+str(knodo1)+'''a'''+str(knodo2)+'''?
				>> '''))
				otroRango = int(input('''
					¿Quieres seleccionar otro rango?
					1 - Si
					* - No
					>> '''))
				if otroRango == 1:
					k[knodo1:knodo2] = valorNuevoK
				else:
					k[knodo1:knodo2] = valorNuevoK
					break
		else:
			k = np.random.rand(N+2)

	else:
		valorK = int(input('''
			Selecciona el único valor de k para toda la barra:
			>> '''))
		k = np.ones(N+2)*valorK

	print('k = ',k)
	plt.plot(k, 'o')
	plt.title('Distribución del valor de k sobre los nodos')
	plt.show()
	return k


def sol_Dirichlet(tolerancia, a, b, N, Tmax, boundA, boundB, k):
	h = (b-a)/(N+1)
	dt = 0.01
	Nt = int(Tmax / dt)
	x = np.linspace(a,b,N+2)
	u = np.zeros(N+2)
	u[0] = boundA
	u[N+1] = boundB
	f = np.copy(u[1:N+1])
	uold = np.copy(u)

	# Construccion de la matriz
	r = dt / (h*h) 
	A = Laplaciano_Dirichlet(N,r,k)

	sumat = 0.0
	error = []

	fig, (ax1, ax2) = plt.subplots(2,1)
	ax1.plot(x, u, '--k',label='Inicial')
	line, = ax1.plot(x, u)
	label = ax1.text(0.75, -.75, 'Time = {:>8.5f}'.format(0), ha='center', va='center', fontsize=12)
	# Ciclo en el tiempo, desde 0 hasta Nt-1
	def animation_frame(n):
	#for n in range(Nt+1):
		time_step = n * dt
		f = np.copy(u[1:N+1])
		uold = np.copy(u)
		sumat = 0.0
		# Solucion en el espacio con Euler hacia atras: IMPLICITO
		t1_start = time.perf_counter()    
		f[0] += r * k[1] * boundA
		f[N-1] += r * k[N+1] * boundB
		u[1:N+1] = np.linalg.solve(A,f) # Sol. del sistema lineal
		t1_stop = time.perf_counter()
		sumat += (t1_stop - t1_start)

		err = np.sqrt(h) * np.linalg.norm(uold-u)
		error.append(err)
		print("n = ", n, ' Error = %12.10g' % err)
		line.set_ydata(u)
		ax2.plot(error, color='blue')

		t1_start = time.perf_counter()    
		f = np.copy(u[1:N+1])
		uold = np.copy(u)
		t1_stop = time.perf_counter()
		sumat += (t1_stop - t1_start)

		if (err < tolerancia):
			plt.pause(0)
		label.set_text('Error = {:>8.5e}\n Time = {:>8.5f}'.format(err, time_step))
		ax1.set_title('Euler Implicito')
	animation = FuncAnimation(fig, func = animation_frame, frames = Nt+1, interval = 0.001, repeat=False)

	plt.xlabel('$x$')
	plt.ylabel('$u(x)$')
	plt.grid()
	plt.legend()
	plt.show()

#sol_Dirichlet(tolerancia, a, b, N, Tmax, boundA, boundB, k)
def sol_Neumman(tolerancia, a, b, N, Tmax, boundA, boundB, k):
	h = (b-a)/(N+1)
	dt = 0.01
	Nt = int(Tmax / dt)
	x = np.linspace(a,b,N+2)
	u = np.zeros(N+2)
	u[0] = boundA
	u[N+1] = boundB+boundA
	f = np.copy(u[1:N+1])
	uold = np.copy(u)

	# Construccion de la matriz
	r = dt / (h*h) 
	A = Laplaciano_Neumman(N,r,k)

	sumat = 0.0
	error = []

	fig, (ax1, ax2) = plt.subplots(2,1)
	ax1.plot(x, u, '--k',label='Inicial')
	line, = ax1.plot(x, u)
	label = ax1.text(0.75, -.75, 'Time = {:>8.5f}'.format(0), ha='center', va='center', fontsize=12)
	# Ciclo en el tiempo, desde 0 hasta Nt-1
	def animation_frame(n):
	#for n in range(Nt+1):
		time_step = n * dt
		f = np.copy(u[1:N+1])
		uold = np.copy(u)
		sumat = 0.0
		# Solucion en el espacio con Euler hacia atras: IMPLICITO
		t1_start = time.perf_counter()    
		f[0] += r * k[1] * boundA
		f[N-1] += r * k[N+1] * boundB * h
		u[1:N+1] = np.linalg.solve(A,f) # Sol. del sistema lineal
		t1_stop = time.perf_counter()
		sumat += (t1_stop - t1_start)

		err = np.sqrt(h) * np.linalg.norm(uold-u)
		error.append(err)
		print("n = ", n, ' Error = %12.10g' % err)
		line.set_ydata(u)
		ax2.plot(error, color='blue')

		t1_start = time.perf_counter()    
		f = np.copy(u[1:N+1])
		uold = np.copy(u)
		t1_stop = time.perf_counter()
		sumat += (t1_stop - t1_start)

		if (err < tolerancia):
			plt.pause(0)
		label.set_text('Error = {:>8.5e}\n Time = {:>8.5f}'.format(err, time_step))
		ax1.set_title('Euler Implicito')
	animation = FuncAnimation(fig, func = animation_frame, frames = Nt+1, interval = 0.001, repeat=False)
	
	plt.xlabel('$x$')
	plt.ylabel('$u(x)$')
	plt.grid()
	plt.legend()
	plt.show()


#inicializando el programa:
tipoFrontera = tipo_Fronteras()
tolerancia, a, b, N, Tmax, boundA, boundB = lectura_Variables(tipoFrontera)
k = lectura_K(N)
if tipoFrontera == 'Dirichlet':
	sol_Dirichlet(tolerancia, a, b, N, Tmax, boundA, boundB, k)
elif tipoFrontera == 'Neumman':
	sol_Neumman(tolerancia, a, b, N, Tmax, boundA, boundB, k)