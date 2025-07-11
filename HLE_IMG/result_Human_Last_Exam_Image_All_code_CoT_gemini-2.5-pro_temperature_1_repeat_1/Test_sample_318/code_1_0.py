import numpy as np
from scipy.integrate import solve_ivp

def solve_epidemiological_enigma():
    """
    This function solves the problem by simulating the disease model to find threshold times,
    using a pre-determined solution to the parameter puzzle, and calculating the final cipher.
    """

    # --- Part 1: Epidemiological Threshold Time Calculation ---

    # Define model parameters
    params = {
        'beta_h': 0.1, 'a_i': 6.0, 'a_r': 21.0, 'f_s': 0.2, 'mu': 1.0/45625.0,
        'mu_n': 1.0/2100.0, 'mu_s': 1.0/600.0, 'mu_h': 1.0/2400.0,
        'c_h': 1.0, 'c_l': 1.0, 'r_b': 0.0
    }
    quarantine_start = 60
    quarantine_end = 116
    beta_n_normal = 0.5
    beta_s_normal = 0.5
    beta_n_quarantine = 0.125
    beta_s_quarantine = 0.125
    
    # Define initial conditions for the state variables
    y0 = [
        999999.0,  # S(0)
        0.0,       # E(0)
        1.0,       # In(0)
        0.0,       # Is(0)
        0.0,       # H(0)
        0.0,       # R(0)
        0.0,       # D(0)
        2900000.0, # B(0)
        0.0,       # Ch(0)
        0.0        # Cl(0)
    ]

    def ode_system(t, y, p):
        """The system of differential equations for the disease model."""
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        # Adjust contact rates based on quarantine period
        if quarantine_start <= t <= quarantine_end:
            beta_n, beta_s = beta_n_quarantine, beta_s_quarantine
        else:
            beta_n, beta_s = beta_n_normal, beta_s_normal
            
        # Total population for infection denominator, per problem statement
        T = max(1.0, E + In + Is + R + S)
        Is_unhosp = max(0, Is - H)
        
        infection_force = (beta_s * S * Is_unhosp / T) + (p['beta_h'] * H * S / T) + (beta_n * S * In / T)

        # Equations
        dSdt = -infection_force - p['mu'] * S
        dEdt = infection_force - E * (1/p['a_i'] + p['mu'])
        from_E_to_I = E / p['a_i']
        dIndt = from_E_to_I * (1 - p['f_s']) - In / p['a_r'] - p['mu_n'] * In
        dIsdt = from_E_to_I * p['f_s'] - Is / p['a_r'] - p['mu_h'] * H - p['mu_s'] * Is_unhosp
        
        if H < B:
            new_hospitalizations = min(max(0, B - H), from_E_to_I * p['f_s'])
        else:
            new_hospitalizations = 0
            
        dHdt = new_hospitalizations - H / p['a_r'] - p['mu_h'] * H
        dRdt = (In + Is) / p['a_r'] - p['mu'] * R
        dDdt_val = p['mu_h'] * H + p['mu_s'] * Is_unhosp + p['mu_n'] * In
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt_val + In + Is)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt_val, dBdt, dChdt, dCldt]

    # Define events to find t_n1 and t_n2
    def event_n1(t, y, p): return y[1] - y[2]  # E - In = 0
    event_n1.direction = 1  # Trigger when crossing from negative to positive
    
    def event_n2(t, y, p): return y[8] - y[6]  # Ch - D = 0
    event_n2.direction = 1  # Trigger when crossing from negative to positive

    # Solve the ODE system
    t_span = [0, 365]
    sol = solve_ivp(
        lambda t, y: ode_system(t, y, params),
        t_span, y0, method='RK45', dense_output=True, events=(event_n1, event_n2)
    )

    # Extract event times and convert from days to hours
    t_n1_hours = (sol.t_events[0][0] * 24.0) if sol.t_events[0].size > 0 else -1
    
    t_n2_event_times_days = [t for t in sol.t_events[1] if t > 1e-6]
    t_n2_hours = (t_n2_event_times_days[0] * 24.0) if t_n2_event_times_days else -1

    # --- Part 2 & 3: Parameter Puzzle Solution and Final Calculation ---

    # Solution to the "Tracing Paper Parameter Puzzle"
    p_map = {
        1: 6, 2: 5, 3: 9, 4: 14, 5: 7, 6: 3, 7: 12, 8: 2, 9: 8
    }
    
    # Calculate X0
    X0 = sum(n * p_n for n, p_n in p_map.items())

    # Calculate the final answer
    final_answer = t_n2_hours * (X0 - t_n1_hours)

    # Output the components of the final equation as requested
    print(f"{t_n2_hours} * ({X0} - {t_n1_hours}) = {final_answer}")
    
    # Output the final answer in the specified format
    print(f"<<<{final_answer}>>>")

solve_epidemiological_enigma()