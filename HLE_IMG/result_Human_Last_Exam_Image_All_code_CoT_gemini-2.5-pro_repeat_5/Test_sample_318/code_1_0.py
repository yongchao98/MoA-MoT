import numpy as np
from scipy.integrate import solve_ivp

def solve_epidemiological_enigma():
    """
    This function solves the entire problem by simulating the model,
    determining the threshold times, and calculating the final cipher.
    """

    # --- Part 1: Model Simulation and Threshold Time Calculation ---

    # Model Parameters (Nominal)
    params = {
        'a_i': 6.0, 'a_r': 21.0, 'f_s': 0.2, 'mu': 1.0/45625.0,
        'mu_n': 1.0/2100.0, 'mu_s': 1.0/600.0, 'mu_h': 1.0/2400.0,
        'beta_h': 0.1, 'c_h': 1.0, 'c_l': 1.0, 'B0': 2900000.0,
        'beta_norm_val': 0.5, 'beta_quar_val': 0.125,
        'quar_start': 60, 'quar_end': 116,
    }

    # Initial Conditions
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    # S, E, I_n, I_s, H, R, D, C_h, C_l

    def seir_model(t, y, p):
        S, E, I_n, I_s, H, R, D, C_h, C_l = [max(0, val) for val in y]

        if p['quar_start'] <= t <= p['quar_end']:
            beta_n = p['beta_quar_val']
            beta_s = p['beta_quar_val']
        else:
            beta_n = p['beta_norm_val']
            beta_s = p['beta_norm_val']

        T = S + E + I_n + I_s + R
        if T == 0:
            return np.zeros(9)

        I_s = max(I_s, H)

        infection_from_In = beta_n * S * I_n / T
        infection_from_Is_unhosp = beta_s * S * (I_s - H) / T
        infection_from_H = p['beta_h'] * S * H / T
        total_infection_rate = infection_from_In + infection_from_Is_unhosp + infection_from_H

        dSdt = -total_infection_rate - p['mu'] * S
        dEdt = total_infection_rate - E * (1/p['a_i'] + p['mu'])
        
        new_severe_rate = E * p['f_s'] / p['a_i']
        new_normal_rate = E * (1 - p['f_s']) / p['a_i']

        dIndt = new_normal_rate - I_n / p['a_r'] - p['mu_n'] * I_n

        h_in = 0
        if H < p['B0']:
            # Assuming B-H is a rate in people/day
            h_in = min(p['B0'] - H, new_severe_rate)
        
        dHdt = h_in - H / p['a_r'] - p['mu_h'] * H

        death_rate_severe = p['mu_h'] * H + p['mu_s'] * (I_s - H)
        dIsdt = new_severe_rate - I_s / p['a_r'] - death_rate_severe
        
        dRdt = (I_n + I_s) / p['a_r'] - p['mu'] * R

        death_rate_normal = p['mu_n'] * I_n
        dDdt = death_rate_severe + death_rate_normal
        
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + I_n + I_s)

        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dChdt, dCldt]

    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], 365 * 24 + 1)
    sol = solve_ivp(seir_model, t_span, y0, args=(params,), dense_output=True, t_eval=t_eval, method='RK45')

    S, E, I_n, I_s, H, R, D, C_h, C_l = sol.y

    # Find t_n1: smallest t > 0 where E(t) > I_n(t)
    diff_E_In = E - I_n
    cross_indices_n1 = np.where(np.diff(np.sign(diff_E_In)) > 0)[0]
    idx1 = cross_indices_n1[0]
    t_cross_n1_days = np.interp(0, [diff_E_In[idx1], diff_E_In[idx1+1]], [sol.t[idx1], sol.t[idx1+1]])
    t_n1 = t_cross_n1_days * 24

    # Find t_n2: smallest t > 0 where C_h(t) > D(t)
    diff_Ch_D = C_h - D
    cross_indices_n2 = np.where(np.diff(np.sign(diff_Ch_D)) > 0)[0]
    idx2 = cross_indices_n2[0]
    t_cross_n2_days = np.interp(0, [diff_Ch_D[idx2], diff_Ch_D[idx2+1]], [sol.t[idx2], sol.t[idx2+1]])
    t_n2 = t_cross_n2_days * 24

    # --- Part 2: Parameter Puzzle Solution ---
    # Based on the deductive analysis explained in the thinking process.
    p = {
        1: 6,  # Plot 1: S vs f_s
        2: 9,  # Plot 2: S vs beta_h
        3: 1,  # Plot 3: S vs mu
        4: 5,  # Plot 4: S vs a_i
        5: 7,  # Plot 5: C_l vs c_l
        6: 2,  # Plot 6: S vs mu_s
        7: 12, # Plot 7: C_h vs c_h
        8: 8,  # Plot 8: S vs mu_h
        9: 3   # Plot 9: R vs mu_n
    }

    # --- Part 3: Final Cipher Calculation ---
    X_0 = sum(n * p_n for n, p_n in p.items())

    final_answer = t_n2 * (X_0 - t_n1)

    print(f"Threshold time t_n1 (E > I_n): {t_n1:.4f} hours")
    print(f"Threshold time t_n2 (C_h > D): {t_n2:.4f} hours")
    print("Deduced parameter mapping {plot: parameter_id}:")
    print(p)
    print(f"Calculated X_0 = sum(n * p_n): {X_0}")
    print("\nFinal calculation:")
    print(f"t_n2 * (X_0 - t_n1) = {t_n2:.4f} * ({X_0} - {t_n1:.4f})")
    print(f"= {t_n2:.4f} * ({X_0 - t_n1:.4f})")
    print(f"= {final_answer:.4f}")
    print("\nFinal Answer:")
    print(f"<<<{final_answer:.1f}>>>")

solve_epidemiological_enigma()