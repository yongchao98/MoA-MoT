import numpy as np
from scipy.integrate import odeint
import math

def solve_epidemiological_enigma():
    """
    Solves the multi-part epidemiological modeling problem.
    """

    # Part 1: Simulate the ODE model to find t_n1 and t_n2

    # --- Model Parameters ---
    # Using a dictionary for clarity
    params = {
        'a_i': 6.0, 'a_r': 21.0, 'f_s': 0.2, 'mu': 1/45625.0,
        'mu_n': 1/2100.0, 'mu_s': 1/600.0, 'mu_h': 1/2400.0,
        'beta_h': 0.1, 'c_h': 1.0, 'c_l': 1.0, 'r_b': 0.0,
        'beta_n_normal': 0.5, 'beta_s_normal': 0.5,
        'beta_n_q': 0.125, 'beta_s_q': 0.125,
        'q_start': 60.0, 'q_end': 116.0
    }

    # --- ODE System Definition ---
    def model(y, t, p):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        # Time-dependent contact rates for quarantine period
        if p['q_start'] <= t <= p['q_end']:
            beta_n = p['beta_n_q']
            beta_s = p['beta_s_q']
        else:
            beta_n = p['beta_n_normal']
            beta_s = p['beta_s_normal']

        # Total population for force of infection calculation
        T = max(1e-9, S + E + In + Is + R)
        
        # Non-hospitalized severely symptomatic individuals
        Is_minus_H = max(0, Is - H)
        
        # Total force of infection
        infections_term = (beta_s * S * Is_minus_H / T) + \
                          (p['beta_h'] * H * S / T) + \
                          (beta_n * S * In / T)

        # Equations of the model
        dSdt = -infections_term - p['mu'] * S
        dEdt = infections_term - E * (1/p['a_i'] + p['mu'])
        dIndt = E * (1 - p['f_s']) / p['a_i'] - In / p['a_r'] - p['mu_n'] * In
        dIsdt = (E * p['f_s'] / p['a_i']) - (Is / p['a_r']) - (p['mu_h'] * H) - (p['mu_s'] * Is_minus_H)
        
        if H < B:
          influx_H = min(B - H, E * p['f_s'] / p['a_i'])
        else:
          influx_H = 0
        influx_H = max(0, influx_H)
        dHdt = influx_H - H / p['a_r'] - p['mu_h'] * H
        
        dRdt = (In + Is) / p['a_r'] - p['mu'] * R
        dDdt = p['mu_h'] * H + p['mu_s'] * Is_minus_H + p['mu_n'] * In
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + In + Is)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # --- Initial Conditions ---
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]

    # --- Simulation Time ---
    t_end = 365
    t = np.linspace(0, t_end, t_end * 24 + 1) # Hourly resolution

    # --- Solve ODE System ---
    sol = odeint(model, y0, t, args=(params,))
    S, E, In, Is, H, R, D, B, Ch, Cl = sol.T

    # --- Calculate t_n1 and t_n2 ---
    try:
        idx_n1 = np.where(E > In)[0][0]
        t_day_n1 = t[idx_n1]
        t_n1 = math.ceil(t_day_n1 * 24)
    except IndexError:
        t_n1 = -1 # Should not happen in this model

    try:
        idx_n2 = np.where(Ch > D)[0][0]
        t_day_n2 = t[idx_n2]
        t_n2 = math.ceil(t_day_n2 * 24)
    except IndexError:
        t_n2 = -1 # Should not happen in this model
        
    # Part 2: Solve the puzzle to find X_0
    # Based on qualitative analysis of parameter effects on the system dynamics.
    # The varied parameter for each plot n is p_n.
    p = {
        1: 9,   # Plot 1: S vs beta_h (strong effect on final size)
        2: 2,   # Plot 2: S vs mu_s (large effect on final size)
        3: 1,   # Plot 3: S vs mu (small effect, curves are close)
        4: 3,   # Plot 4: S vs mu_n (medium effect on final size)
        5: 12,  # Plot 5: Ch vs c_h (S-shaped scaling plot)
        6: 8,   # Plot 6: S vs mu_h (smaller effect than mu_n)
        7: 6,   # Plot 7: C_l vs f_s (S-shaped with crossover)
        8: 11,  # Plot 8: S vs r_b (assuming non-zero r_b variation)
        9: 5,   # Plot 9: H vs a_i (bell-shaped with horizontal shift)
    }

    X_0 = sum(n * p[n] for n in range(1, 10))

    # Part 3: Calculate the final result
    final_answer = t_n2 * (X_0 - t_n1)

    # Print out the components of the calculation as required.
    print(f"Calculated t_n1 (first hour when E > I_n): {t_n1}")
    print(f"Calculated t_n2 (first hour when C_h > D): {t_n2}")
    print("Identified parameter mapping {p_n}:")
    for i in range(1, 10):
        print(f"p_{i}: {p[i]}")
    print(f"Calculated X_0 = sum(n * p_n): {X_0}")
    print(f"Final Calculation: t_n2 * (X_0 - t_n1)")
    print(f"= {t_n2} * ({X_0} - {t_n1})")
    print(f"= {t_n2} * ({X_0 - t_n1})")
    print(f"= {final_answer}")
    return final_answer

if __name__ == '__main__':
    result = solve_epidemiological_enigma()
    print(f"\nFinal Answer: {result}")
    print(f"<<<{result}>>>")
