import numpy as np # type: ignore
from scipy.integrate import solve_ivp # type: ignore # type: ignore
import matplotlib.pyplot as plt

#Definomos las constantes

m2=2.0
G=1
lambda_m=0.1
m1_0=2.0

#Definiimos las funciones

def omega(t):
    return 2.0 * np.sin(t)

def m1(t):
    return m1_0 * np.exp(-lambda_m * t)

def mu(t):
    return (m1(t) * m2) / (m1(t) + m2)

def V(t,r):
    return -G * m1(t) * m2 /( r**2)

#Definimos una funcion numerica para calcular derivadas

def derivative(f, x, dx=1e-6, *args):
    return (f(x + dx, *args) - f(x - dx, *args)) / (2 * dx)

#Sacamos las derivadas parciales necesarias

def dV_r(t,r):
    return derivative(lambda r_: V(t, r_), r)

def dmu_t(t):
    return derivative(mu, t)

#Definimos el sistema de ecuaciones

def equations(t, y):
    """Sistema de ODE:"""
    r, v = y
    mu_i = mu(t)
    drdt = v
    dvdt = r * omega(t)**2 - dV_r(t, r)/mu_i - (dmu_t(t)/mu_i) * V(t, r)

    return [drdt, dvdt]

#Condiciones iniciales: r(0) = r0, r'(0) = v0
r0 = 1.0
v0 = 0.0
initial_conditions = [r0, v0]

# Puntos en el tiempo y tiempo total
t_span = (0, 10)
t_eval = np.linspace(0, 10, 1000)

# Solucion del ode
solution = solve_ivp(equations, t_span, initial_conditions, t_eval=t_eval)

# Resultados
r_sol = solution.y[0] 
v_sol = solution.y[1]

# Save time (t) and radial position (r) to a text file
data = np.column_stack((solution.t, solution.y[0]))
np.savetxt('radial_motion.txt', data, header='time\t radial_position')



    