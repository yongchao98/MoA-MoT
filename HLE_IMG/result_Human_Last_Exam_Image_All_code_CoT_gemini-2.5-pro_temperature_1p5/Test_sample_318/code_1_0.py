import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_epidemic_enigma():
    """
    Solves the multi-part epidemiological problem.
    1. Simulates the ODE model to find t_n1 and t_n2.
    2. Defines the parameter mapping p_n found through qualitative analysis.
    3. Calculates X0 and the final answer.
    """

    # --- Part 1: Epidemiological Threshold Time Calculation ---

    # Model Parameters
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

    # Initial conditions
    S0 = 999999.0
    I_n0 = 1.0
    B0 = 2900000.0
    E0, I_s0, H0, R0, D0, C_h0, C_l0 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    y0 = [S0, E0, I_n0, I_s0, H0, R0, D0, B0, C_h0, C_l0]

    def get_betas(t):
        if 60 <= t <= 116:
            return 0.125, 0.125
        else:
            return 0.5, 0.5

    def model_as_written(t, y):
        S, E, I_n, I_s, H, R, D, B, C_h, C_l = y
        beta_n_t, beta_s_t = get_betas(t)
        
        T = max(0, E + I_n + I_s + R + S)
        if T == 0:
            return np.zeros(10)

        contact_s = beta_s_t * S * max(0, I_s - H) / T
        contact_h = beta_h * H * S / T
        contact_n = beta_n_t * S * I_n / T
        
        dSdt = -contact_s - contact_h - contact_n - mu * S
        dEdt = contact_s + contact_h + contact_n - E * (1/a_i + mu)
        dI_ndt = E * (1 - f_s) / a_i - I_n / a_r - mu_n * I_n
        dI_sdt = E * f_s / a_i - I_s / a_r - mu_h * H - mu_s * (I_s - H)
        
        hospitalization_inflow = 0
        if H < B:
            hospitalization_inflow = min(B - H, E * f_s / a_i)
        dHdt = hospitalization_inflow - H / a_r - mu_h * H

        dRdt = (I_n + I_s) / a_r - mu * R
        dDdt = mu_h * H + mu_s * (I_s - H) + mu_n * I_n
        dBdt = r_b * B
        dC_hdt = c_h * H
        dC_ldt = c_l * (dDdt + I_n + I_s) 

        return [dSdt, dEdt, dI_ndt, dI_sdt, dHdt, dRdt, dDdt, dBdt, dC_hdt, dC_ldt]

    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], t_span[1] * 24 * 10 + 1) # High resolution
    sol = solve_ivp(model_as_written, t_span, y0, method='RK45', t_eval=t_eval)

    # Find t_n1: smallest t > 0 (in hours) such that E > I_n
    t_n1 = -1
    for i in range(1, len(sol.t)):
        if sol.y[1, i] > sol.y[2, i]:
            t_n1 = sol.t[i] * 24
            break

    # Find t_n2: smallest t > 0 (in hours) such that C_h > D
    t_n2 = -1
    for i in range(1, len(sol.t)):
        if sol.y[8, i] > sol.y[6, i]:
            t_n2 = sol.t[i] * 24
            break
            
    # --- Part 2 & 3: Cipher Calculation ---

    # Parameter mapping {p1, p2, ..., p9} based on analysis
    p = {
        1: 5,   # a_i
        2: 14,  # q_l
        3: 15,  # q_f
        4: 13,  # q_s
        5: 12,  # c_h
        6: 6,   # f_s
        7: 8,   # mu_h
        8: 2,   # mu_s
        9: 9    # beta_h
    }

    # Calculate X0
    X0 = sum(n * p[n] for n in range(1, 10))

    # Calculate final answer
    final_answer = t_n2 * (X0 - t_n1)

    print(f"t_n1 (first time E > I_n) = {t_n1:.4f} hours")
    print(f"t_n2 (first time C_h > D) = {t_n2:.4f} hours")
    print("Parameter mapping:")
    print("p = {1:5, 2:14, 3:15, 4:13, 5:12, 6:6, 7:8, 8:2, 9:9}")
    print(f"X0 = 1*5 + 2*14 + 3*15 + 4*13 + 5*12 + 6*6 + 7*8 + 8*2 + 9*9 = {X0}")
    print(f"Final calculation: {t_n2:.4f} * ({X0} - {t_n1:.4f}) = {final_answer:.4f}")
    
    # Final answer in the required format
    print(f"<<<{final_answer}>>>")


solve_epidemic_enigma()