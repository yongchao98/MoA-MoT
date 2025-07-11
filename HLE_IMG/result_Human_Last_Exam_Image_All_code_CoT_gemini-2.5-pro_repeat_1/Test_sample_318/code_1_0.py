import numpy as np
from scipy.integrate import solve_ivp
import math

# Plan:
# 1. Define the parameters and the system of ordinary differential equations (ODEs) as specified in the problem.
# 2. Solve the ODE system numerically to get the time evolution of all compartments.
# 3. From the solution, find the specific time points t_n1 and t_n2.
#    - t_n1: The first time (in hours) when the Exposed population (E) exceeds the Normally Infected population (I_n).
#    - t_n2: The first time (in hours) when the cumulative healthcare cost (C_h) exceeds the cumulative deceased population (D).
# 4. Define the vector of parameter identifiers {p_1, ..., p_9} based on a qualitative analysis of the provided plots. This part involves deducing which parameter variation corresponds to each plot. This deduction is done beforehand and the result is hardcoded.
# 5. Calculate X_0 = sum(n * p_n for n=1 to 9).
# 6. Calculate the final answer using the formula: t_n2 * (X_0 - t_n1).
# 7. Print all the components of the final formula (t_n1, t_n2, X_0) and the final result.

def solve_epidemiological_enigma():
    """
    Solves the entire problem by implementing the model, finding the thresholds,
    using the pre-determined parameter mappings, and calculating the final cipher.
    """

    # --- Part 1: Epidemiological Threshold Time Calculation ---

    # Step 1: Define parameters and ODE system
    # Model parameters
    beta_h = 0.1
    a_i = 6.0
    a_r = 21.0
    f_s = 0.2
    mu = 1.0 / 45625.0
    mu_n = 1.0 / 2100.0
    mu_s = 1.0 / 600.0
    mu_h = 1.0 / 2400.0
    c_h = 1.0
    c_l = 1.0
    r_b = 0.0

    # Initial conditions
    y0 = [
        999999.0,  # S(0)
        0.0,       # E(0)
        1.0,       # I_n(0)
        0.0,       # I_s(0)
        0.0,       # H(0)
        0.0,       # R(0)
        0.0,       # D(0)
        2900000.0, # B(0)
        0.0,       # C_h(0)
        0.0        # C_l(0)
    ]

    def get_betas(t):
        """Returns time-dependent contact rates."""
        if 60.0 <= t <= 116.0:
            return 0.125, 0.125  # beta_n, beta_s during quarantine
        else:
            return 0.5, 0.5  # beta_n, beta_s outside quarantine

    def model(t, y):
        """The system of ODEs."""
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        T = S + E + In + Is + R
        if T <= 0: T = 1.0

        beta_n, beta_s = get_betas(t)
        
        Is_minus_H_non_neg = np.max([0.0, Is - H])
        Is_minus_H = Is - H
        
        # Transmission dynamics
        infection_rate = (beta_s * S * Is_minus_H_non_neg / T) + \
                         (beta_h * H * S / T) + \
                         (beta_n * S * In / T)
        
        # ODEs as per problem description
        dSdt = -infection_rate - mu * S
        dEdt = infection_rate - E * (1.0/a_i + mu)
        dIndt = E * (1.0 - f_s) / a_i - In / a_r - mu_n * In
        
        new_severe_cases_rate = E * f_s / a_i
        
        inflow_H = 0.0
        if H < B:
            inflow_H = np.min([B - H, new_severe_cases_rate])
            inflow_H = np.max([0.0, inflow_H])
            
        dHdt = inflow_H - H / a_r - mu_h * H
        dIsdt = new_severe_cases_rate - Is / a_r - mu_h * H - mu_s * Is_minus_H
        dRdt = (In + Is) / a_r - mu * R
        dDdt = mu_h * H + mu_s * Is_minus_H + mu_n * In
        dBdt = r_b * B
        dChdt = c_h * H
        dCldt = c_l * (dDdt + In + Is)

        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Step 2: Solve the ODE system
    t_span = [0, 365]
    sol = solve_ivp(model, t_span, y0, dense_output=True, method='RK45', rtol=1e-6, atol=1e-6)

    # Step 3: Find t_n1 and t_n2
    # Create a fine time grid to find the crossing points accurately
    t_fine = np.linspace(t_span[0], t_span[1], (t_span[1] - t_span[0]) * 100)
    y_fine = sol.sol(t_fine)

    E_vals, In_vals = y_fine[1], y_fine[2]
    Ch_vals, D_vals = y_fine[8], y_fine[9]

    # Find t_n1: smallest time t > 0 where E(t) > I_n(t)
    try:
        idx_n1 = np.where(E_vals > In_vals)[0][0]
        # Linear interpolation for a more precise crossing time
        t1, t2 = t_fine[idx_n1 - 1], t_fine[idx_n1]
        e1, in1 = E_vals[idx_n1 - 1], In_vals[idx_n1 - 1]
        e2, in2 = E_vals[idx_n1], In_vals[idx_n1]
        t_n1_days = t1 + (in1 - e1) * (t2 - t1) / (e2 - e1 - in2 + in1)
        t_n1 = t_n1_days * 24.0  # Convert to hours
    except IndexError:
        t_n1 = -1 # Should not happen based on problem dynamics

    # Find t_n2: smallest time t > 0 where C_h(t) > D(t)
    try:
        # D grows faster initially, so we look for the first crossing after t=0
        idx_n2 = np.where(Ch_vals > D_vals)[0][0]
        # Linear interpolation
        t1, t2 = t_fine[idx_n2 - 1], t_fine[idx_n2]
        ch1, d1 = Ch_vals[idx_n2 - 1], D_vals[idx_n2 - 1]
        ch2, d2 = Ch_vals[idx_n2], D_vals[idx_n2]
        t_n2_days = t1 + (d1 - ch1) * (t2 - t1) / (ch2 - ch1 - d2 + d1)
        t_n2 = t_n2_days * 24.0 # Convert to hours
    except IndexError:
        t_n2 = -1 # Should not happen

    # --- Part 2: Tracing Paper Parameter Puzzle ---

    # Step 4: Define the parameter identifiers {p_n}
    # These values are deduced from qualitative analysis of the plots and the model equations.
    # p_n is the identifier for the parameter varied in plot n.
    # The mapping is: μ:1, μ_s:2, μ_n:3, a_i:5, f_s:6, c_l:7, μ_h:8, β_h:9, r_b:11, c_h:12, q_s:13, q_l:14, q_f:15
    p = {
        1: 14,  # Plot 1: S vs q_l (quarantine length)
        2: 5,   # Plot 2: S vs a_i (incubation period)
        3: 1,   # Plot 3: S vs mu (baseline mortality)
        4: 13,  # Plot 4: S vs q_s (quarantine start)
        5: 9,   # Plot 5: D vs beta_h (hospital contact rate)
        6: 3,   # Plot 6: S vs mu_n (normal infected mortality)
        7: 7,   # Plot 7: C_l vs c_l (lost productivity cost)
        8: 8,   # Plot 8: S vs mu_h (hospitalized mortality)
        9: 6    # Plot 9: H vs f_s (fraction severe)
    }

    # Step 5: Calculate X_0
    X_0 = sum(n * p[n] for n in range(1, 10))

    # --- Part 3: The Epidemiological Enigma's Cipher ---
    
    # Step 6: Calculate the final answer
    final_answer = t_n2 * (X_0 - t_n1)

    # Step 7: Print the components and the final result
    print(f"t_n1 (hours) = {t_n1}")
    print(f"t_n2 (hours) = {t_n2}")
    print(f"X_0 = {X_0}")
    print(f"Final equation: {t_n2:.4f} * ({X_0} - {t_n1:.4f})")
    print(f"Final Answer = {final_answer}")
    print(f"<<<{final_answer:.1f}>>>")


solve_epidemiological_enigma()