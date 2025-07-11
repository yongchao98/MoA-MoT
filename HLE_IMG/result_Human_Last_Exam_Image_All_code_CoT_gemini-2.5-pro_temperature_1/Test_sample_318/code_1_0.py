import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_epidemiological_enigma():
    """
    This function solves the multi-part epidemiological modeling problem.
    It first calculates the threshold times t_n1 and t_n2 by numerically solving the ODE system.
    Then, it uses a pre-determined mapping of plots to parameters to calculate X_0.
    Finally, it computes the final answer using the given formula.
    """

    # --- Part 1: Numerical Simulation and Threshold Time Calculation ---

    # Define the model parameters based on the problem description
    params = {
        'beta_h': 0.1,
        'a_i': 6.0,
        'a_r': 21.0,
        'f_s': 0.2,
        'mu': 1/45625,
        'mu_n': 1/2100,
        'mu_s': 1/600,
        'mu_h': 1/2400,
        'c_h': 1.0,
        'c_l': 1.0,
        'r_b': 0.0,
        'beta_n_base': 0.5,
        'beta_s_base': 0.5,
        'q_start': 60,
        'q_end': 116,
        'q_factor': 0.125 / 0.5,  # Quarantine reduction factor
    }

    # Define the initial conditions for the state variables
    # State vector: [S, E, I_n, I_s, H, R, D, B, C_h, C_l]
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

    # Define the system of Ordinary Differential Equations (ODEs)
    def model(t, y, p):
        S, E, In, Is, H, R, D, B, Ch, Cl = y

        # Determine contact rates, which are time-dependent due to quarantine
        if p['q_start'] <= t <= p['q_end']:
            beta_n = p['beta_n_base'] * p['q_factor']
            beta_s = p['beta_s_base'] * p['q_factor']
        else:
            beta_n = p['beta_n_base']
            beta_s = p['beta_s_base']
        
        # Total population T(t) for transmission denominator
        T = max(1e-9, S + E + In + Is + R)

        # Transmission terms
        S_to_E = (beta_n * In + beta_s * max(0, Is - H) + p['beta_h'] * H) * S / T
        
        # Progression from Exposed to Infected
        E_to_I = E / p['a_i']

        # Derivatives for each compartment
        dSdt = -S_to_E - p['mu'] * S
        dEdt = S_to_E - E * (1/p['a_i'] + p['mu'])
        dIndt = E_to_I * (1 - p['f_s']) - In / p['a_r'] - p['mu_n'] * In
        
        # Hospitalization inflow is limited by bed availability
        new_hospitalizations = 0
        if H < B:
            new_hospitalizations = min(B - H, E_to_I * p['f_s'])
        
        dHdt = new_hospitalizations - H / p['a_r'] - p['mu_h'] * H
        dIsdt = E_to_I * p['f_s'] - Is / p['a_r'] - p['mu_h'] * H - p['mu_s'] * max(0, Is - H)
        dRdt = (In + Is) / p['a_r'] - p['mu'] * R
        dDdt = p['mu_h'] * H + p['mu_s'] * max(0, Is - H) + p['mu_n'] * In
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + In + Is)

        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Define event functions to find t_n1 and t_n2
    def event_E_gt_In(t, y, p): return y[1] - y[2] # E - I_n = 0
    event_E_gt_In.direction = 1  # Find crossing from negative to positive

    def event_Ch_gt_D(t, y, p): return y[8] - y[6] # C_h - D = 0
    event_Ch_gt_D.direction = 1  # Find crossing from negative to positive

    # Solve the ODE system
    t_span = (0, 365)
    sol = solve_ivp(
        fun=model,
        t_span=t_span,
        y0=y0,
        args=(params,),
        events=(event_E_gt_In, event_Ch_gt_D),
        dense_output=True,
        method='Radau'
    )

    # Extract event times and convert from days to hours
    t_n1 = sol.t_events[0][0] * 24 if len(sol.t_events[0]) > 0 else -1
    t_n2 = sol.t_events[1][0] * 24 if len(sol.t_events[1]) > 0 else -1

    # --- Part 2: Tracing Paper Puzzle Solution ---
    
    # The unique parameter varied in each plot {p_1, ..., p_9} was determined by
    # analyzing the qualitative effects of each parameter on the model dynamics.
    # The mapping from plot number to parameter ID is given below.
    # Parameter IDs: mu:1, mu_s:2, mu_n:3, a_i:5, f_s:6, c_l:7, mu_h:8, beta_h:9, r_b:11, c_h:12, q_s:13, q_l:14, q_f:15
    p_n_map = {
        1: 13,  # Plot 1: Varies q_s (quarantine start time)
        2: 5,   # Plot 2: Varies a_i (incubation period)
        3: 1,   # Plot 3: Varies mu (baseline mortality)
        4: 6,   # Plot 4: Varies f_s (fraction severe)
        5: 12,  # Plot 5: Varies c_h (healthcare cost per bed)
        6: 8,   # Plot 6: Varies mu_h (hospitalized mortality)
        7: 2,   # Plot 7: Varies mu_s (severe mortality)
        8: 3,   # Plot 8: Varies mu_n (normal mortality)
        9: 9    # Plot 9: Varies beta_h (hospitalized contact rate)
    }

    # --- Final Calculation ---

    # Calculate X_0 = sum_{n=1 to 9} n * p_n
    X_0 = sum(n * p_n_map[n] for n in range(1, 10))

    # Calculate the final answer: t_n2 * (X_0 - t_n1)
    final_answer = t_n2 * (X_0 - t_n1)

    # Print the components of the final equation as requested
    print("Final Equation: t_n2 * (X_0 - t_n1)")
    print(f"t_n1 (time when E > I_n, in hours) = {t_n1}")
    print(f"t_n2 (time when C_h > D, in hours) = {t_n2}")
    print(f"X_0 (sum of n * p_n) = {X_0}")
    print(f"Result = {t_n2} * ({X_0} - {t_n1})")
    print(f"Final Answer = {final_answer}")

# Execute the main function
solve_epidemiological_enigma()