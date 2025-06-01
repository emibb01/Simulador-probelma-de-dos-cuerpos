def run_two_body_simulation(m1, m2, simulation_time=1000):
    """Simulacion principal del problema de dos cuerpos"""

    import numpy as np # type: ignore
    from scipy.integrate import solve_ivp # type: ignore # type: ignore
    import matplotlib.pyplot as plt

    #Definomos las constantes

    G=1 # constante gravitacional normalizada
    d = 1.0 #separacion inicial

    def distance(r1, r2, theta1, theta2):
        return np.sqrt(r1**2 + r2**2 - 2*r1*r2*np.cos(theta1 - theta2)) #Definimos la distancia entre cuerpos

    def V(r1,r2,theta1,theta2):
        return -G * m1 * m2 /distance(r1,r2,theta1, theta2) #Potencial gravitatorio

    #Definimos una funcion numerica para calcular derivadas

    def derivative(f, x, dx=1e-6, *args):
        return (f(x + dx, *args) - f(x - dx, *args)) / (2 * dx) #Definimos una funcion para las derivadas numericas

    def dV_r1(r1,r2,theta1,theta2):
        return derivative(lambda r1_: V(r1_, r2,theta1,theta2), r1) #Obtenemos la derivada del potencial respecto a r1

    def dV_r2(r1,r2,theta1,theta2):
        return derivative(lambda r2_: V(r1, r2_,theta1,theta2), r2) #Obtenemos la derivada del respecto a r2


    #Definimos el sistema de ecuaciones

    def equations(t, y):
        """Sistema de ODE del problema de dos cuerpos:"""
        r1, v_r1, theta1, omega1, r2, v_r2, theta2, omega2 = y
        """Primer cuerpo"""
        dr1dt = v_r1
        dv_r1dt = r1 * omega1**2 - (1/m1) * dV_r1(r1, r2, theta1, theta2)
        dtheta1dt = omega1
        domega1dt = -2 * v_r1 * omega1 / r1  # Aceleracion angular
        """Sistema del segundo cuerpo"""
        dr2dt = v_r2
        dv_r2dt = r2 * omega2**2 - (1/m2) * dV_r2(r1, r2, theta1, theta2)
        dtheta2dt = omega2
        domega2dt = -2 * v_r2 * omega2 / r2  # Aceleracion angular
        return [dr1dt, dv_r1dt, dtheta1dt, domega1dt, 
                dr2dt, dv_r2dt, dtheta2dt, domega2dt]

    

    # Velocidad orbital para orbitaws circulares, esta es la velocidad inicial supuesta
    v1 = np.sqrt(G * m2**2 / (d * (m1 + m2))) 
    v2 = np.sqrt(G * m1**2 / (d * (m1 + m2)))

    # Postion inicial de los cuerpos
    x1, y1 = -d, 0.0
    x2, y2 = d, 0.0

    # Velocidad inicial de ambos cuerpos
    vx1, vy1 = 0.0, v1
    vx2, vy2 = 0.0, -v2

    # Convertimos las coordenadas polares a cartesianas
    r1_0 = np.sqrt(x1**2 + y1**2)
    theta1_0 = np.arctan2(y1, x1)
    vr1_0 = (x1 * vx1 + y1 * vy1) / r1_0
    omega1_0 = (x1 * vy1 - y1 * vx1) / r1_0**2

    r2_0 = np.sqrt(x2**2 + y2**2)
    theta2_0 = np.arctan2(y2, x2)
    vr2_0 = (x2 * vx2 + y2 * vy2) / r2_0
    omega2_0 = (x2 * vy2 - y2 * vx2) / r2_0**2

    # Vector de condiciones iniciales [r1(0), vr1(0), theta1(0), omega1(0), r2(0), vr2(0), theta2(0), omega2(0)]
    y0 = [r1_0, vr1_0, theta1_0, omega1_0,
        r2_0, vr2_0, theta2_0, omega2_0]
    
    # Definimos el tiempo de la simulacion
    t_span = (0, simulation_time)
    t_eval = np.linspace(*t_span, 1000)

    # Resolvemos ODE
    solution = solve_ivp(equations, t_span, y0, t_eval=t_eval, method='DOP853', rtol=1e-8)

    # Extramos la posicion de los cuerpos de la solucion
    r1 = solution.y[0]
    theta1 = solution.y[2]
    r2 = solution.y[4]
    theta2 = solution.y[6]

    # Convertimos las posiciones polares a cartesianas
    x1 = r1 * np.cos(theta1)
    y1 = r1 * np.sin(theta1)
    x2 = r2 * np.cos(theta2)
    y2 = r2 * np.sin(theta2)

    # Grafica de las orbitas
    plt.figure(figsize=(8, 8))
    plt.plot(x1, y1, label='Cuerpo 1')
    plt.plot(x2, y2, label='Cuerpo 2')
    plt.scatter(0, 0, color='k', marker='+', label='COM')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Orbitas de dos cuerpos')
    plt.legend()
    plt.grid()
    plt.axis('equal')
    plt.show()

import tkinter as tk
from tkinter import ttk

def run_simulation_from_gui():
    try:
        m1 = float(mass1_entry.get())
        m2 = float(mass2_entry.get())
        time = float(time_entry.get())
        
        # Close the interface window
        root.destroy()
        
        # Run the simulation
        run_two_body_simulation(m1, m2, time)
        
    except ValueError:
        error_label.config(text="Please enter valid numbers!")

# Create the interface
root = tk.Tk()
root.title("Masa de los cuerpos a simular")

# Make the window resizable and set minimum size
root.minsize(300, 200)
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

# Main frame
main_frame = ttk.Frame(root, padding="20")
main_frame.grid(row=0, column=0, sticky="nsew")
main_frame.columnconfigure(1, weight=1)

# Mass 1 input
ttk.Label(main_frame, text="Masa 1 (kg):").grid(row=0, column=0, sticky="w", pady=5)
mass1_entry = ttk.Entry(main_frame)
mass1_entry.grid(row=0, column=1, sticky="ew", pady=5)
mass1_entry.insert(0, "150.0")

# Mass 2 input
ttk.Label(main_frame, text="Masa 2 (kg):").grid(row=1, column=0, sticky="w", pady=5)
mass2_entry = ttk.Entry(main_frame)
mass2_entry.grid(row=1, column=1, sticky="ew", pady=5)
mass2_entry.insert(0, "200.0")

# Simulation time input
ttk.Label(main_frame, text="Tiempo de simulacion:").grid(row=2, column=0, sticky="w", pady=5)
time_entry = ttk.Entry(main_frame)
time_entry.grid(row=2, column=1, sticky="ew", pady=5)
time_entry.insert(0, "10.0")

# Run button
run_button = ttk.Button(main_frame, text="Correr simulacion", command=run_simulation_from_gui)
run_button.grid(row=3, column=0, columnspan=2, pady=15)

# Error label
error_label = ttk.Label(main_frame, text="", foreground="red")
error_label.grid(row=4, column=0, columnspan=2)

root.mainloop()