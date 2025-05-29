import numpy as np

class Planeta(object):
    global G, M, m
    G = 1
    M = 1
    m = 1

    def __init__(self, condicion_inicial, alpha=0):
        self.y_actual = condicion_inicial
        self.t_actual = 0
        self.alpha = alpha

    def ecuacion_de_movimiento(self):
        x, y, vx, vy = self.y_actual
        fx = lambda x, y, t: (2 * self.alpha * G * M * x) / ((x**2 + y**2)**2) - (G * M * x) / ((np.sqrt(x**2 + y**2))**3)
        fy = lambda x, y, t: (2 * self.alpha * G * M * y) / ((x**2 + y**2)**2) - (G * M * y) / ((np.sqrt(x**2 + y**2))**3)
        return [vx, vy, fx, fy]

    def avanza_verlet(self, dt):
        t0 = self.t_actual
        x0, y0, vx0, vy0 = self.y_actual
        fx = self.ecuacion_de_movimiento()[2]
        fy = self.ecuacion_de_movimiento()[3]
        xn = x0 + vx0 * dt + (fx(x0, y0, t0) * (dt**2)) / 2.0
        yn = y0 + vy0 * dt + (fy(x0, y0, t0) * (dt**2)) / 2.0
        vxn = vx0 + ((fx(x0, y0, t0) + fx(xn, yn, t0 + dt)) * dt) / 2.0
        vyn = vy0 + ((fy(x0, y0, t0) + fy(xn, yn, t0 + dt)) * dt) / 2.0
        self.y_actual = xn, yn, vxn, vyn

    def energia_total(self):
        x0, y0, vx0, vy0 = self.y_actual
        E = 0.5 * m * (vx0**2 + vy0**2) + (self.alpha * G * M * m) / (x0**2 + y0**2) - (G * M * m) / (np.sqrt(x0**2 + y0**2))
        return E

def main_solution(condicion_inicial, alpha, dt, n):
    condicion_inicial = np.array(condicion_inicial)
    dt = float(dt)
    n = int(n)
    
    P = Planeta(condicion_inicial, alpha)
    
    for _ in range(n):
        P.avanza_verlet(dt)
    
    final_energy = P.energia_total()
    
    return final_energy

# Further adjusted parameters
condicion_inicial = [1.0, 0.0, 0.0, 0.96]
alpha = 0.22
dt = 0.01
n = 1000

final_energy = main_solution(condicion_inicial, alpha, dt, n)
print(final_energy)