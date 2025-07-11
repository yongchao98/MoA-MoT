import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_and_calculate():
    """
    This function implements the full solution logic.
    It first solves the ODE system to find t_n1 and t_n2.
    Then, it uses the pre-determined solution to the parameter puzzle
    to calculate X_0. Finally, it computes and prints the required cipher.
    """
    # PART 1: Find t_n1 and t_n2 by solving the ODE system
    # -----------------------------------------------------

    # Model parameters and initial conditions from the problem description
    S0 = 999999
    In0 = 1
    B0 = 2900000
    y0 = [S0, 0, In0, 0, 0, 0, 0, B0, 0, 0] # S, E, In, Is, H, R, D, B, Ch, Cl

    params_list = [
        0.1,         # beta_h: contact rate for hospitalized
        6.0,         # a_i: incubation period
        21.0,        # a_r: infectious period
        0.2,         # f_s: fraction of infected developing severe symptoms
        1/45625.0,   # mu: baseline mortality rate
        1/2100.0,    # mu_n: mortality rate for normally infected
        1/600.0,     # mu_s: mortality rate for severely infected
        1/2400.0,    # mu_h: mortality rate for hospitalized
        1.0,         # c_h: healthcare cost unit
        1.0,         # c_l: lost productivity cost unit
        0.0          # r_b: hospital beds change rate
    ]
    
    # Time-dependent contact rates with quarantine period
    def beta_n_s(t):
        return 0.125 if 60 <= t <= 116 else 0.5

    # Definition of the system of ODEs
    def model(t, y, params):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        beta_h, a_i, a_r, f_s, mu, mu_n, mu_s, mu_h, c_h, c_l, r_b = params

        # Total living population as defined in the model
        T = max(0, S + E + In + Is + R)
        T = T if T > 0 else 1

        beta_t = beta_n_s(t)
        Is_minus_H = max(0, Is - H) # Non-hospitalized severe cases
        
        # Combined infection rate from all infectious groups
        infection_term = (beta_t * S * Is_minus_H + beta_h * H * S + beta_t * S * In) / T
        
        # Derivatives for each compartment based on the problem's equations
        dSdt = -infection_term - mu * S
        dEdt = infection_term - E * (1/a_i + mu)
        
        new_severe_flux = E * f_s / a_i
        
        dIndt = E * (1 - f_s) / a_i - In / a_r - mu_n * In
        
        hosp_in = min(B - H, new_severe_flux) if H < B else 0
        dHdt = hosp_in - H / a_r - mu_h * H
        
        dIsdt = new_severe_flux - Is / a_r - mu_h * H - mu_s * Is_minus_H
        dDdt = mu_h * H + mu_s * Is_minus_H + mu_n * In
        dRdt = (In + Is) / a_r - mu * R
        dBdt = r_b * B
        dChdt = c_h * H
        dCldt = c_l * (dDdt + In + Is)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Solve the ODE system
    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], 365 * 24 + 1) # Hourly evaluation for precision
    sol = solve_ivp(model, t_span, y0, args=(params_list,), dense_output=True, method='RK45', rtol=1e-8, atol=1e-8)
    y_eval = sol.sol(t_eval)

    # Extract solution arrays for relevant compartments
    E_sol, In_sol, D_sol, Ch_sol = y_eval[1], y_eval[2], y_eval[6], y_eval[8]
    
    # Find t_n1: smallest time t (in hours) where E > In
    idx_n1 = np.where(E_sol > In_sol)[0][0]
    t_n1 = t_eval[idx_n1] * 24

    # Find t_n2: smallest time t (in hours) where Ch > D
    idx_n2 = np.where(Ch_sol > D_sol)[0][0]
    t_n2 = t_eval[idx_n2] * 24

    # PART 2: Define X_0 based on the solution to the parameter puzzle
    # ----------------------------------------------------------------
    # The mapping from plot number 'n' to parameter ID 'p_n' is hardcoded here.
    # It is the result of qualitative analysis of the model's behavior under parameter variations.
    p = {
        1: 9,   # Plot 1 -> beta_h
        2: 2,   # Plot 2 -> mu_s
        3: 8,   # Plot 3 -> mu_h
        4: 3,   # Plot 4 -> mu_n
        5: 12,  # Plot 5 -> c_h
        6: 6,   # Plot 6 -> f_s
        7: 15,  # Plot 7 -> q_f (quarantine factor)
        8: 1,   # Plot 8 -> mu
        9: 5    # Plot 9 -> a_i
    }
    X_0 = sum(n * p_n for n, p_n in p.items())

    # PART 3: Calculate and print the final answer
    # ---------------------------------------------
    final_answer = t_n2 * (X_0 - t_n1)

    # Print the final equation with all its numerical components, as requested.
    print(f"{t_n2} * ({X_0} - {t_n1}) = {final_answer}")

# Execute the function to get the final answer.
solve_and_calculate()