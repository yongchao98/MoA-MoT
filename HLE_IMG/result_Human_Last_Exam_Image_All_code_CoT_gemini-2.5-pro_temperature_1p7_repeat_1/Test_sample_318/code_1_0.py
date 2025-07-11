import numpy as np
from scipy.integrate import solve_ivp
import math

def main():
    """
    Main function to solve the entire problem.
    """

    # --- Step 1: Define Model and Parameters ---

    # Nominal model parameters
    params = {
        'beta_h': 0.1, 'beta_n_normal': 0.5, 'beta_n_quarantine': 0.125,
        'beta_s_normal': 0.5, 'beta_s_quarantine': 0.125,
        'a_i': 6.0, 'a_r': 21.0, 'f_s': 0.2, 'mu': 1.0 / 45625.0,
        'mu_n': 1.0 / 2100.0, 'mu_s': 1.0 / 600.0, 'mu_h': 1.0 / 2400.0,
        'c_h': 1.0, 'c_l': 1.0, 'r_b': 0.0, 'q_start': 60, 'q_end': 116,
    }

    # Initial conditions for the state vector y = [S, E, In, Is, H, R, D, B, Ch, Cl]
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]

    def seir_model(t, y, p):
        """
        Defines the system of differential equations.
        """
        S, E, In, Is, H, R, D, B, Ch, Cl = y

        # Time-dependent contact rates for quarantine period
        if p['q_start'] <= t <= p['q_end']:
            beta_n = p['beta_n_quarantine']
            beta_s = p['beta_s_quarantine']
        else:
            beta_n = p['beta_n_normal']
            beta_s = p['beta_s_normal']

        T = max(0, E + In + Is + R + S)
        if T == 0: T = 1 # Avoid division by zero

        Is_nh = max(0, Is - H) # Non-hospitalized severely symptomatic
        lambda_val = (beta_s * Is_nh + p['beta_h'] * H + beta_n * In) / T

        dSdt = -lambda_val * S - p['mu'] * S
        dEdt = lambda_val * S - (1/p['a_i'] + p['mu']) * E
        dIndt = (1 - p['f_s']) * E / p['a_i'] - In / p['a_r'] - p['mu_n'] * In
        
        # Rate of new hospitalizations, limited by bed availability
        if H < B:
            h_inflow = min(max(0, B - H), (p['f_s'] * E / p['a_i']))
        else:
            h_inflow = 0
            
        dIsdt = p['f_s'] * E / p['a_i'] - Is / p['a_r'] - p['mu_h'] * H - p['mu_s'] * Is_nh
        dHdt = h_inflow - H / p['a_r'] - p['mu_h'] * H
        dRdt = (In + Is) / p['a_r'] - p['mu'] * R
        
        dDdt = p['mu_h'] * H + p['mu_s'] * Is_nh + p['mu_n'] * In
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + In + Is)

        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # --- Step 2: Calculate t_n1 and t_n2 ---

    # Define event functions to find the precise time thresholds
    event_E_gt_In = lambda t, y, p: y[1] - y[2]  # E - In
    event_Ch_gt_D = lambda t, y, p: y[8] - y[6]  # Ch - D

    # Solve the ODE system with event detection
    sol = solve_ivp(seir_model, (0, 365), y0, args=(params,), dense_output=True,
                    events=[event_E_gt_In, event_Ch_gt_D], method='RK45')

    # Extract event times and convert from days to hours, then apply ceiling
    t_n1 = -1
    if sol.t_events[0].size > 0:
        t_n1 = math.ceil(sol.t_events[0][0] * 24)

    t_n2 = -1
    if sol.t_events[1].size > 0:
        t_n2 = math.ceil(sol.t_events[1][0] * 24)

    # --- Step 3: Use deduced parameter mappings for p_n ---
    
    # Based on qualitative analysis of the plot shapes and parameter effects.
    # Parameter IDs: mu:1, mu_s:2, mu_n:3, a_i:5, f_s:6, c_l:7, mu_h:8, 
    #                beta_h:9, c_h:12, q_s:13
    p_map = {
        1: 9,   # Plot 1: S(t) vs beta_h (amplifying, kink at same time)
        2: 6,   # Plot 2: S(t) vs f_s (dampening, strong effect)
        3: 1,   # Plot 3: S(t) vs mu (dampening, very weak effect)
        4: 13,  # Plot 4: S(t) vs q_s (amplifying, kink time varies)
        5: 12,  # Plot 5: C_h(t) vs c_h (increasing, scaled version)
        6: 3,   # Plot 6: S(t) vs mu_n (dampening, strong effect)
        7: 7,   # Plot 7: C_l(t) vs c_l (increasing, scaled version)
        8: 2,   # Plot 8: S(t) vs mu_s (dampening, moderate effect)
        9: 5    # Plot 9: E(t) vs a_i (humped shape, peak time shifts)
    }

    # --- Step 4: Calculate X0 ---

    X0 = sum(n * p_map[n] for n in range(1, 10))

    # --- Step 5: Final Calculation ---

    final_result = t_n2 * (X0 - t_n1)

    # Print the equation with all the components as requested
    print(f"Final Calculation: t_n2 * (X0 - t_n1)")
    print(f"= {t_n2} * ({X0} - {t_n1})")
    print(f"= {final_result}")
    
    # Final answer in the required format
    print(f"<<<{final_result}>>>")

if __name__ == "__main__":
    main()