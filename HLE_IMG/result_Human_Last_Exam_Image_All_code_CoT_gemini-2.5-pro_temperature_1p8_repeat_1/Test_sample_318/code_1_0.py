import numpy as np
from scipy.integrate import solve_ivp

def solve_epidemic_model():
    """
    Solves the epidemic model to find threshold times, deciphers the parameter puzzle,
    and computes the final answer.
    """

    # Part 1: Numerical Simulation and Threshold Calculation
    def seird_model(t, y, params):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        a_i = params['a_i']
        a_r = params['a_r']
        f_s = params['f_s']
        mu = params['mu']
        mu_n = params['mu_n']
        mu_s = params['mu_s']
        mu_h = params['mu_h']
        beta_h = params['beta_h']
        beta_n_nom = params['beta_n_nom']
        beta_s_nom = params['beta_s_nom']
        r_b = params['r_b']
        c_h = params['c_h']
        c_l = params['c_l']

        if 60 <= t <= 116:
            beta_n = beta_n_nom / 4.0
            beta_s = beta_s_nom / 4.0
        else:
            beta_n = beta_n_nom
            beta_s = beta_s_nom

        T = max(0, S + E + In + Is + R)
        if T == 0: T = 1

        Is_non_H = max(0, Is - H)
        infection_force = (beta_s * S * Is_non_H / T) + (beta_h * H * S / T) + (beta_n * S * In / T)

        dS_dt = -infection_force - mu * S
        dE_dt = infection_force - E * (1/a_i + mu)
        dIn_dt = E * (1-f_s)/a_i - In/a_r - mu_n * In
        
        hosp_inflow = 0
        if H < B:
            hosp_inflow = min(B - H, E * f_s / a_i)
        
        dIs_dt = E * f_s / a_i - Is / a_r - mu_h * H - mu_s * (Is - H)
        dH_dt = hosp_inflow - H / a_r - mu_h * H
        dR_dt = (In + Is) / a_r - mu * R
        
        dD_dt_val = mu_h * H + mu_s * (Is - H) + mu_n * In
        dCl_dt_val = c_l * (dD_dt_val + In + Is)
        dCh_dt_val = c_h * H
        dB_dt_val = r_b * B
        
        return [dS_dt, dE_dt, dIn_dt, dIs_dt, dH_dt, dR_dt, dD_dt_val, dB_dt_val, dCh_dt_val, dCl_dt_val]

    params_nom = {
        'a_i': 6.0, 'a_r': 21.0, 'f_s': 0.2, 'mu': 1/45625, 'mu_n': 1/2100,
        'mu_s': 1/600, 'mu_h': 1/2400, 'beta_h': 0.1, 'beta_n_nom': 0.5,
        'beta_s_nom': 0.5, 'r_b': 0, 'c_h': 1, 'c_l': 1
    }
    y0 = [999999, 0, 1, 0, 0, 0, 0, 2900000, 0, 0]
    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], 365*24 + 1)
    
    sol_nom = solve_ivp(seird_model, t_span, y0, args=(params_nom,), dense_output=True, t_eval=t_eval)
    
    t_sol = sol_nom.t
    E, In, D, Ch = sol_nom.y[1], sol_nom.y[2], sol_nom.y[6], sol_nom.y[8]

    # Find t_n1: smallest t > 0 where E > In
    idx_n1 = np.where(E > In)[0][0]
    t_n1_hours = t_sol[idx_n1] * 24

    # Find t_n2: smallest t > 0 where Ch > D
    idx_n2 = np.where(Ch[1:] > D[1:])[0][0]
    t_n2_hours = t_sol[idx_n2 + 1] * 24

    # Part 2: Tracing Paper Puzzle Solution
    p = {
        1: 9,  # beta_h
        2: 6,  # f_s
        3: 8,  # mu_h
        4: 3,  # mu_n
        5: 12, # c_h
        6: 2,  # mu_s
        7: 7,  # c_l
        8: 1,  # mu
        9: 5   # a_i
    }

    # Part 3: Final Calculation
    X0 = sum(n * p[n] for n in range(1, 10))
    
    final_answer = t_n2_hours * (X0 - t_n1_hours)

    print("--- Epidemiological Enigma's Cipher ---")
    print(f"Equation: t_n2 * (X_0 - t_n1)")
    print(f"Value for t_n1 (first time in hours E > I_n): {t_n1_hours:.1f}")
    print(f"Value for t_n2 (first time in hours C_h > D): {t_n2_hours:.1f}")
    print(f"Value for X_0 (from puzzle): {X0}")
    print(f"Calculation: {t_n2_hours:.1f} * ({X0} - {t_n1_hours:.1f})")
    print(f"Final Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")


solve_epidemic_model()