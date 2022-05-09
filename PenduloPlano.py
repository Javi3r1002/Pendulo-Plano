"""
Universidad del Valle de Guatemala           4/5/2022
Mecánica I                Javier Mejía Alecio (20304)

                Parcial II
Con este programa se busca determinar el ángulo de un péndulo plano a través del tiempo utilizando el método de Euler y el modelo
de un artículo. Además se comparan experimentalmente
"""
#Se importan las librerías a utilizar
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import mpmath as mpm
from sklearn import metrics
import math
from math import pi
#Se definene las variables a utilizar
ai = pi/6 # ángulo inicial
s = ai
l = 0.97 
g = 9.8 
Delta = 0.05 # Ya que la frecuencia de medición son de 20hz

#Se exportan las mediciones de Capstone
Data = pd.read_csv('Pendulo_Plano.csv')

Data = Data.dropna()

print(Data.columns.values)

#Se utilizan los datos de la corrida en la que se utilizó un ángulo de 30°
Data = Data[['Time (s) Run #8', 'Angle (rad) Run #8']]
Data = Data[Data["Time (s) Run #8"] > 7.2]
Tiempo = Data.iloc[:,0]-Data.iloc[0,0]
# Se crea una array de numpy para tener los mismos tiempos en los modelos y en los datos experimentales
Tiempo = Tiempo.to_numpy()
T = np.arange(Tiempo[0],Tiempo[-1],Delta)
N = len(T)
for i in range(N):
    T[i] = T[i]-Tiempo[0]


#Función Jacobiana

Yj=[]
a = (math.sin(ai/2))
K = mpm.ellipk(a**2)
w = np.sqrt(g/l)

def E(t):
    G = mpm.ellipfun('sn', K-w*t,a**2)
    return(G)


for t in T:
    y = 2*math.asin(a*(E(t)))
    Yj.append(y)

# Método de Euler
Ve = []
W = 0
k = s
for i in T:
    Ai = k
    W += ((-g/l))*math.sin(Ai)*Delta
    k += W*Delta
    Ve.append(k)

Ve = np.array(Ve)
Yj = np.array(Yj)
residuales1 = np.array(Data.iloc[:,1])-Yj
residuales2 = np.array(Data.iloc[:,1])-Ve
plt.boxplot(residuales1);
plt.title("Diagrama de caja y bigotes de residuos entre datos experimentales y el modelo del artículo:")
plt.ylabel('Ángulo')
plt.show()

plt.boxplot(residuales2);
plt.title("Diagrama de caja y bigotes de residuos entre datos experimentales y el método de Euler:")
plt.ylabel('Ángulo')
plt.show()

plt.plot(T,Ve,'b.-',T,Yj,'r-', T, Data.iloc[:,1],'g')
plt.legend(['Euler','Artículo', 'Experimental'])
#plt.plot(T, Ve, 'b.-', T, Yj, 'r-')
#plt.legend(['Euler','Artículo'])
#plt.title("Gráficos ángulo vs tiempo") # Titulo del gráfico
plt.xlabel('t') #Título eje x
plt.ylabel('θ') #Título eje y
plt.show()