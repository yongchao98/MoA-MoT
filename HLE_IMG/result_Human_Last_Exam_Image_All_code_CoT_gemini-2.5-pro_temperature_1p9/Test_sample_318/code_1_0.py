import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_puzzle():
    """
    This function solves the entire problem by:
    1. Defining and solving the epidemiological model to find t_n1 and t_n2.
    2. Using a pre-determined mapping of plots to parameters to calculate X0.
    3. Computing the final expression t_n2 * (X0 - t_n1).
    """

    # --- Part 1: Model Definition ---

    # Define model parameters
    params = {
        'a_i': 6, 'a_r': 21, 'f_s': 0.2, 'mu': 1/45625, 'mu_n': 1/2100,
        'mu_s': 1/600, 'mu_h': 1/2400, 'beta_h': 0.1, 'r_b': 0, 'c_h': 1, 'c_l': 1
    }

    # Define time-dependent contact rates
    def beta_n(t): return 0.125 if 60 <= t <= 116 else 0.5
    def beta_s(t): return 0.125 if 60 <= t <= 116 else 0.5

    # Define initial conditions for [S, E, I_n, I_s, H, R, D, B, C_h, C_l]
    y0 = [999999., 0., 1., 0., 0., 0., 0., 2900000., 0., 0.]

    # Define the system of differential equations
    def epidemic_model(t, y, p):
        S, E, I_n, I_s, H, R, D, B, C_h, C_l = y
        
        T = max(1e-9, E + I_n + I_s + R + S)
        I_s_unhosp = max(0, I_s - H)
        
        b_n_t = beta_n(t)
        b_s_t = beta_s(t)
        
        dSdt = - (b_s_t * S * I_s_unhosp) / T - (p['beta_h'] * H * S) / T - (b_n_t * S * I_n) / T - p['mu'] * S
        dEdt = (b_s_t * S * I_s_unhosp) / T + (p['beta_h'] * H * S) / T + (b_n_t * S * I_n) / T - E * (1/p['a_i'] + p['mu'])
        dIndt = E * (1 - p['f_s']) / p['a_i'] - I_n / p['a_r'] - p['mu_n'] * I_n
        
        h_in_potential = E * p['f_s'] / p['a_i']
        h_in = min(max(0, B - H), h_in_potential) if H < B else 0
        
        dIsdt = h_in_potential - I_s / p['a_r'] - p['mu_h'] * H - p['mu_s'] * I_s_unhosp
        dHdt = h_in - H / p['a_r'] - p['mu_h'] * H
        dRdt = (I_n + I_s) / p['a_r'] - p['mu'] * R
        dDdt = p['mu_h'] * H + p['mu_s'] * I_s_unhosp + p['mu_n'] * I_n
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + I_n + I_s)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # --- Part 1 continued: Find Threshold Times t_n1 and t_n2 ---

    # Define event functions for the solver
    def event_tn1(t, y, p): return y[1] - y[2] # E > I_n
    event_tn1.direction = 1

    def event_tn2(t, y, p): return y[8] - y[6] # C_h > D
    event_tn2.direction = 1

    # Solve the ODE system
    sol = solve_ivp(
        fun=epidemic_model, t_span=[0, 365], y0=y0, args=(params,),
        events=(event_tn1, event_tn2), dense_output=True, rtol=1e-9, atol=1e-9
    )

    # Extract event times in days and convert to hours
    t_n1 = sol.t_events[0][0] * 24
    t_n2 = sol.t_events[1][0] * 24

    # --- Part 2: Define Parameter Mappings ---
    # Based on qualitative analysis matching simulated plot characteristics to the puzzle image.
    p_identifiers = {
        1: 5,  # Plot 1: S(t) vs a_i (delay effect)
        2: 9,  # Plot 2: S(t) vs beta_h (high impact)
        3: 2,  # Plot 3: S(t) vs mu_s (large severity effect)
        4: 3,  # Plot 4: S(t) vs mu_n (moderate severity effect)
        5: 12, # Plot 5: C_h(t) vs c_h (perfect scaling)
        6: 6,  # Plot 6: S(t) vs f_s (strong global effect)
        7: 7,  # Plot 7: C_l(t) vs c_l (perfect scaling)
        8: 1,  # Plot 8: S(t) vs mu (minor effect)
        9: 8,  # Plot 9: D(t) vs mu_h (non-scaling dynamic effect on deaths)
    }

    # --- Part 3: Final Calculation ---

    # Calculate X0
    X0 = sum(n * p_identifiers[n] for n in range(1, 10))

    # Calculate the final result
    final_answer = t_n2 * (X0 - t_n1)

    # Print the equation with substituted values
    print("Final equation with calculated values:")
    print(f"{t_n2} * ({X0} - {t_n1})")
    
    # Print the final result
    print("\nCalculated result:")
    print(final_answer)
    
    return final_answer

final_result = solve_puzzle()
# The final result should be presented inside <<<>>>
# For example: <<< -1089330.14 >>>
# Let's format the output this way for the final submission.
# The calculation in the block will provide the number.
