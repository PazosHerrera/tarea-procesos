#Se importan las librerias necesarias
import numpy as np #Para operaciones con vectores
import math # para funciones matematicas
import random as rnd #generacion de numeros aleatorios
import matplotlib.pyplot as plt #para graficar


#condiciones iniciales
t0=0.0
X0=0.0
Y0=500.0

#parametros
alpha = 2
beta = 1.
gamma = 1.
delta = 1.
Omega= 1.
Steps = 50000
#Se define la Matriz estoquiometrica vista en clase
S =[[1.,0.,-1.,0.],[0.,1.,0.,-1.]]
#Se inicializa los vectores solucion y el tiempo
Y = np.zeros([2,Steps+1])
t = np.zeros(Steps+1)
#Se inicializan las condiciones iniciales
Y[0][0]=X0
Y[1][0]=Y0

#funcion que genera numero aleatorio de distribucion exponencial
def dist_exp(a):
	r = rnd.random() #Se genera un numero aleatorio de distribucion uniforme
	return -(1./a)*math.log(r) #Se tranforma en distribucion exponencial

#Elige aleatoriamente la reaccion definida en el vector
# de propensiones ni y una tasa de reaccion global A
def dist_reaction(ni,A):
	r = rnd.random() #Numero aleatorio de distribucion uniforme
	if (r < ni[0]/A): #Se elige reaccion 0
		return 0
	elif (r < (ni[0] + ni[1])/A): #Se elige reaccion 1
		return 1
	elif (r < (ni[0] + ni[1] + ni[2])/A): #Se elige reaccion 2
		return 2
	elif (r <= (ni[0] + ni[1] + ni[2] + ni[3])/A): #Se elige reaccion 3
		return 3
#ciclo del algoritmo
for i in range(Steps):
	#se actualiza vector de propensiones
	ni = np.array([gamma*Omega, alpha*Y[0][i]*(Y[0][i]-1.)*Y[1][i]/Omega**2, \
	beta*Y[0][i], delta*Y[0][i] ])
	#se actualiza tasa de reaccion global
	a = sum(ni)
	#se elige aleatoriamente tiempo de reaccion
	tau = dist_exp(a)
	#se elige la reaccion que se llevara a cabo
	mu = dist_reaction(ni,a)
	#se actualiza el estado del sistema y se avanza el tiempo
	Y[0][i+1] = Y[0][i] + S[0][mu]
	Y[1][i+1] = Y[1][i] + S[1][mu]
	t[i+1] = t[i] + tau

#se grafica el espacio fase
plt.plot(Y[0],Y[1])
#Se grafica la evolucion en el tiempo
plt.figure()
plt.plot(t,Y[0])
#se muestran las graficas
plt.show()
