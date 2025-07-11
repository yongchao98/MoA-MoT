import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_epidemiological_enigma():
    """
    This function solves the entire problem:
    1. Defines and solves the ODE system for the disease outbreak.
    2. Finds the threshold times t_n1 and t_n2.
    3. Calculates X_0 based on the puzzle solution.
    4. Computes and prints the final answer.
    """

    # --- Model Parameters ---
    a_i = 6.0
    a_r = 21.0
    f_s = 0.2
    mu = 1.0 / 45625.0
    mu_n = 1.0 / 2100.0
    mu_s = 1.0 / 600.0
    mu_h = 1.0 / 2400.0
    beta_h = 0.1
    c_h = 1.0
    c_l = 1.0
    r_b = 0.0

    def get_betas(t):
        if 60 <= t <= 116:
            return 0.125  # beta_n and beta_s
        else:
            return 0.5

    def model(t, y):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        beta_val = get_betas(t)
        
        T = S + E + In + Is + R
        if T <= 0: T = 1

        Is_nh = max(0, Is - H)
        
        lambda_val = (beta_val * Is_nh + beta_h * H + beta_val * In) / T
        
        dSdt = -lambda_val * S - mu * S
        dEdt = lambda_val * S - E * (1/a_i + mu)
        dIndt = E * (1 - f_s) / a_i - In / a_r - mu_n * In
        
        dDdt_term_s = mu_s * Is_nh
        dDdt_term_h = mu_h * H
        
        dIsdt = E * f_s / a_i - Is / a_r - dDdt_term_h - dDdt_term_s

        H_increase = 0
        if H < B:
            H_increase = min(B - H, E * f_s / a_i)
        dHdt = H_increase - H / a_r - mu_h * H
        
        dRdt = (In + Is) / a_r - mu * R
        dDdt = dDdt_term_h + dDdt_term_s + mu_n * In
        dBdt = r_b * B
        dChdt = c_h * H
        dCldt = c_l * (dDdt + In + Is)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # --- Initial Conditions and Simulation ---
    y0 = [999999.0, 0, 1.0, 0, 0, 0, 0, 2900000.0, 0, 0]
    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], 365 * 24 + 1) # Hourly evaluation

    sol = solve_ivp(model, t_span, y0, t_eval=t_eval, method='RK45')
    t_days = sol.t
    E_sol, In_sol, D_sol, Ch_sol = sol.y[1], sol.y[2], sol.y[6], sol.y[8]

    # --- Find t_n1: E(t) > In(t) ---
    t_n1_hours = -1
    # Find the index where the condition is first met
    try:
        idx1 = np.where(E_sol > In_sol)[0][0]
        # Linear interpolation for better precision
        t1, t2 = t_days[idx1-1], t_days[idx1]
        v1, v2 = (E_sol - In_sol)[idx1-1], (E_sol - In_sol)[idx1]
        t_n1_days = t1 - v1 * (t2 - t1) / (v2 - v1)
        t_n1_hours = t_n1_days * 24
    except IndexError:
        t_n1_hours = -1 # Should not happen

    # --- Find t_n2: Ch(t) > D(t) ---
    t_n2_hours = -1
    # D starts growing first, so we look for Ch to overtake it.
    try:
        idx2 = np.where(Ch_sol > D_sol)[0][0]
        # Linear interpolation for better precision
        t1, t2 = t_days[idx2-1], t_days[idx2]
        v1, v2 = (Ch_sol - D_sol)[idx2-1], (Ch_sol - D_sol)[idx2]
        t_n2_days = t1 - v1 * (t2 - t1) / (v2 - v1)
        t_n2_hours = t_n2_days * 24
    except IndexError:
        t_n2_hours = -1 # Should not happen

    # --- Calculate X_0 ---
    p = {1: 5, 2: 9, 3: 6, 4: 2, 5: 12, 6: 3, 7: 7, 8: 8, 9: 1}
    X_0 = sum(n * p[n] for n in range(1, 10))

    # --- Final Calculation ---
    final_answer = t_n2_hours * (X_0 - t_n1_hours)

    print(f"Calculation: t_n2 * (X_0 - t_n1)")
    print(f"Values: {t_n2_hours:.4f} * ({X_0} - {t_n1_hours:.4f})")
    print(f"Result: {final_answer:.4f}")
    
    # Final answer in the required format
    print(f"<<<{final_answer}>>>")


solve_epidemiological_enigma()