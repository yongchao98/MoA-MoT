import numpy as np
from scipy.integrate import solve_ivp

def solve_epidemiological_enigma():
    """
    This function solves the problem by first numerically solving the ODE system
    to find t_n1 and t_n2, then calculating X_0 based on a qualitative analysis
    of the parameter plots, and finally computing the final expression.
    """

    # Part 1: Numerical Solution to find t_n1 and t_n2
    
    # Define model parameters based on the problem description
    params = {
        'beta_h': 0.1,
        'a_i': 6.0,
        'a_r': 21.0,
        'f_s': 0.2,
        'mu': 1.0 / 45625.0,
        'mu_n': 1.0 / 2100.0,
        'mu_s': 1.0 / 600.0,
        'mu_h': 1.0 / 2400.0,
        'c_h': 1.0,
        'c_l': 1.0,
        'r_b': 0.0,
    }

    # Define the time-dependent contact rate due to quarantine
    def get_beta_t(t):
        return 0.125 if 60 <= t <= 116 else 0.5

    # Define the system of differential equations
    def model_odes(t, y, p):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        T_pop = max(0, S + E + In + Is + R)
        if T_pop == 0:
            return np.zeros(10)

        beta_t = get_beta_t(t)

        infection_rate = (beta_t * In + beta_t * max(0, Is - H) + p['beta_h'] * H) * S / T_pop
        
        dSdt = -infection_rate - p['mu'] * S
        dEdt = infection_rate - E * (1/p['a_i'] + p['mu'])
        dIndt = E * (1 - p['f_s']) / p['a_i'] - In / p['a_r'] - p['mu_n'] * In
        dIsdt = E * p['f_s'] / p['a_i'] - Is / p['a_r'] - p['mu_h'] * H - p['mu_s'] * (Is - H)
        
        hospitalization_flow = min(B - H, E * p['f_s'] / p['a_i']) if H < B else 0
        dHdt = hospitalization_flow - H / p['a_r'] - p['mu_h'] * H
        
        dRdt = (In + Is) / p['a_r'] - p['mu'] * R
        dDdt = p['mu_h'] * H + p['mu_s'] * (Is - H) + p['mu_n'] * In
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + In + Is)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Set initial conditions
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]
    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], t_span[1] * 24 * 10) # High-resolution time points
    
    # Solve the ODE system
    sol = solve_ivp(lambda t, y: model_odes(t, y, params), t_span, y0, method='RK45', t_eval=t_eval)

    # Extract required compartments
    t_days = sol.t
    E_vals, In_vals, Ch_vals, D_vals = sol.y[1], sol.y[2], sol.y[8], sol.y[6]

    # Find t_n1: smallest time t (in hours) where E > In
    try:
        idx_n1 = np.where(E_vals > In_vals)[0][0]
        # Interpolate between points for precision
        t1, t2 = t_days[idx_n1-1], t_days[idx_n1]
        v1, v2 = (E_vals - In_vals)[idx_n1-1], (E_vals - In_vals)[idx_n1]
        t_n1 = (t1 - v1 * (t2 - t1) / (v2 - v1)) * 24  # Convert to hours
    except IndexError:
        t_n1 = -1

    # Find t_n2: smallest time t (in hours) where Ch > D
    try:
        idx_n2 = np.where(Ch_vals > D_vals)[0][0]
        # Interpolate between points for precision
        t1, t2 = t_days[idx_n2-1], t_days[idx_n2]
        v1, v2 = (Ch_vals - D_vals)[idx_n2-1], (Ch_vals - D_vals)[idx_n2]
        t_n2 = (t1 - v1 * (t2 - t1) / (v2 - v1)) * 24  # Convert to hours
    except IndexError:
        t_n2 = -1

    # Part 2: Tracing Paper Puzzle Solution to find X_0
    # The parameter identifiers p_n for each plot n are deduced qualitatively.
    p_assignments = {
        1: 8,   # Varied: μ_h (mortality in hospital)
        2: 6,   # Varied: f_s (fraction with severe symptoms)
        3: 15,  # Varied: q_f (quarantine factor)
        4: 9,   # Varied: β_h (contact rate for hospitalized)
        5: 12,  # Varied: c_h (healthcare cost unit)
        6: 5,   # Varied: a_i (incubation period)
        7: 7,   # Varied: c_l (lost productivity cost unit)
        8: 1,   # Varied: μ (baseline mortality rate)
        9: 2    # Varied: μ_s (mortality for severe, not hospitalized)
    }
    
    # Calculate X_0 = Σ(n * p_n)
    X_0 = sum(n * p_assignments[n] for n in range(1, 10))

    # Part 3: Final Calculation and Output
    final_answer = t_n2 * (X_0 - t_n1)

    print("--- Epidemiological Threshold Times ---")
    print(f"t_n1 (first time E > I_n) = {t_n1:.4f} hours")
    print(f"t_n2 (first time C_h > D) = {t_n2:.4f} hours")
    
    print("\n--- Puzzle Solution ---")
    print(f"Parameter assignments p_n: {list(p_assignments.values())}")
    print(f"X_0 = sum(n * p_n) = {X_0}")
    
    print("\n--- Final Answer Calculation ---")
    print(f"Result = t_n2 * (X_0 - t_n1)")
    print(f"Result = {t_n2:.4f} * ({X_0} - {t_n1:.4f})")
    print(f"Result = {final_answer:.4f}")
    
    print(f"\n<<<{final_answer}>>>")

solve_epidemiological_enigma()