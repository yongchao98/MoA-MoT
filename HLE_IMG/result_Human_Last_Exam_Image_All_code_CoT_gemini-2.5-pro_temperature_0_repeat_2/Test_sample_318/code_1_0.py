import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_epidemiological_enigma():
    """
    This function solves the entire problem by:
    1. Calculating X_0 based on the puzzle solution.
    2. Simulating the ODE model to find t_n1 and t_n2.
    3. Computing the final answer.
    """
    # Part 2: Tracing Paper Parameter Puzzle Solution
    # p_n is the identifier of the parameter varied in plot n.
    p = {
        1: 3,   # mu_n
        2: 5,   # a_i
        3: 6,   # f_s
        4: 14,  # q_l
        5: 12,  # c_h
        6: 8,   # mu_h
        7: 7,   # c_l
        8: 2,   # mu_s
        9: 11   # r_b
    }
    X_0 = sum(n * p_n for n, p_n in p.items())

    # Part 1: Epidemiological Threshold Time Calculation
    # Model parameters
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
    In0 = 1.0
    B0 = 2900000.0
    E0, Is0, H0, R0, D0, Ch0, Cl0 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    y0 = [S0, E0, In0, Is0, H0, R0, D0, B0, Ch0, Cl0]

    # System of differential equations
    def seir_model(t, y):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        T = max(0, S + E + In + Is + R)
        if T == 0:
            return np.zeros(10)

        if 60 <= t <= 116:
            beta_n, beta_s = 0.125, 0.125
        else:
            beta_n, beta_s = 0.5, 0.5

        unhosp_Is = max(0, Is - H)
        force_of_infection = (beta_n * In + beta_s * unhosp_Is + beta_h * H) * S / T

        dSdt = -force_of_infection - mu * S
        dEdt = force_of_infection - E * (1/a_i + mu)
        
        new_hospitalizations = E * f_s / a_i
        
        dIndt = E * (1 - f_s) / a_i - In / a_r - mu_n * In
        dIsdt = E * f_s / a_i - Is / a_r - mu_h * H - mu_s * unhosp_Is
        dHdt = new_hospitalizations - H / a_r - mu_h * H
        
        dRdt = (In + Is) / a_r - mu * R
        dDdt = mu_h * H + mu_s * unhosp_Is + mu_n * In
        
        dBdt = r_b * B
        dChdt = c_h * H
        dCldt = c_l * (dDdt + In + Is)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Solve the ODE
    t_span = [0, 200]
    sol = solve_ivp(seir_model, t_span, y0, dense_output=True, method='RK45', atol=1e-9, rtol=1e-6)

    # Evaluate solution on a fine time grid to find crossing points
    t_eval_fine = np.linspace(t_span[0], t_span[1], 20001)
    y_eval = sol.sol(t_eval_fine)
    
    # Find t_n1: smallest time t (in hours) such that E > In
    E_vals, In_vals = y_eval[1], y_eval[2]
    t_n1_idx = np.argmax(E_vals > In_vals)
    t_n1 = t_eval_fine[t_n1_idx] * 24

    # Find t_n2: smallest time t (in hours) such that Ch > D
    D_vals, Ch_vals = y_eval[6], y_eval[8]
    t_n2_idx = np.argmax(Ch_vals > D_vals)
    t_n2 = t_eval_fine[t_n2_idx] * 24

    # Final calculation
    result = t_n2 * (X_0 - t_n1)

    # Print the components of the final calculation
    print(f"X_0 = {X_0}")
    print(f"t_n1 (hours) = {t_n1}")
    print(f"t_n2 (hours) = {t_n2}")
    print(f"Final Answer Calculation: {t_n2:.4f} * ({X_0} - {t_n1:.4f}) = {result:.4f}")
    
solve_epidemiological_enigma()