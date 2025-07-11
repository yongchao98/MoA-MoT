import numpy as np
from scipy.integrate import solve_ivp

def solve_problem():
    """
    This function solves the entire problem by first simulating the ODE model
    to find t_n1 and t_n2, then calculating X0 based on the analysis of
    the plots, and finally computing the required expression.
    """

    # --- Part 1: ODE Simulation for t_n1 and t_n2 ---

    # Model parameters from the problem description
    a_i = 6.0
    a_r = 21.0
    f_s = 0.2
    mu = 1.0 / 45625
    mu_n = 1.0 / 2100
    mu_s = 1.0 / 600
    mu_h = 1.0 / 2400
    beta_h = 0.1
    c_h = 1.0
    c_l = 1.0
    r_b = 0.0

    # Define the system of differential equations
    def model(t, y):
        S, E, I_n, I_s, H, R, D, B, C_h, C_l = y
        
        # Time-dependent contact rates based on quarantine period
        beta_val = 0.125 if 60 <= t <= 116 else 0.5
        
        T = max(S + E + I_n + I_s + R, 1e-9)
        unhospitalized_severe = max(0, I_s - H)

        # Equations as specified in the problem
        dSdt = -(beta_val * S * unhospitalized_severe / T) - (beta_h * H * S / T) - (beta_val * S * I_n / T) - mu * S
        dEdt = (beta_val * S * unhospitalized_severe / T) + (beta_h * H * S / T) + (beta_val * S * I_n / T) - E * (1 / a_i + mu)
        dIndt = E * (1 - f_s) / a_i - I_n / a_r - mu_n * I_n
        dIsdt = E * f_s / a_i - I_s / a_r - mu_h * H - mu_s * unhospitalized_severe
        
        new_hospitalizations = min(B - H, E * f_s / a_i) if H < B else 0
        dHdt = new_hospitalizations - H / a_r - mu_h * H
        
        dRdt = (I_n + I_s) / a_r - mu * R
        dDdt = mu_h * H + mu_s * unhospitalized_severe + mu_n * I_n
        dBdt = r_b * B
        dChdt = c_h * H
        dCldt = c_l * (dDdt + I_n + I_s)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Initial conditions
    y0 = [999999, 0, 1, 0, 0, 0, 0, 2900000, 0, 0]
    t_span = [0, 365]

    # Event definitions for t_n1 and t_n2
    def event_E_gt_In(t, y): return y[1] - y[2]  # E - I_n
    event_E_gt_In.direction = 1
    
    def event_Ch_gt_D(t, y): return y[8] - y[6]  # C_h - D
    event_Ch_gt_D.direction = 1

    # Solve the ODE system
    sol = solve_ivp(model, t_span, y0, events=(event_E_gt_In, event_Ch_gt_D), dense_output=True)

    # Extract event times and convert from days to hours
    t_n1 = sol.t_events[0][0] * 24.0
    t_n2 = sol.t_events[1][0] * 24.0

    # --- Part 2: Parameter Identification and X0 Calculation ---

    # Parameter identifiers based on qualitative analysis
    p_identifiers = {
        1: 2, 2: 9, 3: 5, 4: 3, 5: 12,
        6: 8, 7: 7, 8: 1, 9: 6
    }
    
    X0 = sum(n * p_identifiers[n] for n in range(1, 10))

    # --- Part 3: Final Calculation and Output ---
    
    result = t_n2 * (X0 - t_n1)

    print(f"t_n1 (hours, E > I_n): {t_n1}")
    print(f"t_n2 (hours, C_h > D): {t_n2}")
    print(f"p_n mapping: {p_identifiers}")
    print(f"X0 = sum(n * p_n): {X0}")
    print(f"Final calculation: {t_n2:.4f} * ({X0} - {t_n1:.4f}) = {result:.4f}")

    # Return the final numerical answer for submission
    return result

# Run the full process and get the final number
final_answer = solve_problem()
print(f"<<<{final_answer}>>>")