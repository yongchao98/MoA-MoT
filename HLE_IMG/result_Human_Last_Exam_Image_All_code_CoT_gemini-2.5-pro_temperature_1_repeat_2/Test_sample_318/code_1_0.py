import numpy as np
from scipy.integrate import solve_ivp

def solve_epidemiological_enigma():
    """
    Solves the epidemiological model to find t_n1 and t_n2,
    then calculates the final answer based on the pre-calculated X_0.
    """
    
    # Step 1: Define the model and parameters
    # ------------------------------------------

    # The system of differential equations
    def model(t, y, params):
        S, E, I_n, I_s, H, R, D, B, C_h, C_l = y
        
        # Unpack parameters
        a_i = params['a_i']
        a_r = params['a_r']
        f_s = params['f_s']
        mu = params['mu']
        mu_n = params['mu_n']
        mu_s = params['mu_s']
        mu_h = params['mu_h']
        beta_h = params['beta_h']
        c_h = params['c_h']
        c_l = params['c_l']
        r_b = params['r_b']
        
        # Time-dependent contact rates for quarantine period
        if 60 <= t <= 116:
            beta_n = 0.125
            beta_s = 0.125
        else:
            beta_n = 0.5
            beta_s = 0.5

        # Eq 1: Total population for transmission calculation
        T = max(0, S + E + I_n + I_s + R)
        if T == 0: T = 1.0 # Avoid division by zero
            
        # Intermediate terms for infection rates
        term_s_unhosp = beta_s * S * max(0, I_s - H) / T
        term_s_hosp = beta_h * H * S / T
        term_s_normal = beta_n * S * I_n / T
        total_infection_rate = term_s_unhosp + term_s_hosp + term_s_normal

        # System of ODEs
        # Eq 2: S'
        dSdt = -total_infection_rate - mu * S
        # Eq 3: E'
        dEdt = total_infection_rate - E * (1/a_i + mu)
        # Eq 4: I_n'
        dIndt = E * (1 - f_s) / a_i - I_n / a_r - mu_n * I_n
        
        # Eq 6: H'
        if H < B:
            new_hospitalizations = min(B - H, E * f_s / a_i)
        else:
            new_hospitalizations = 0
        dHdt = new_hospitalizations - H / a_r - mu_h * H
        
        # Eq 5: I_s'
        dIsdt = E * f_s / a_i - I_s / a_r - mu_h * H - mu_s * max(0, I_s - H)
        # Eq 7: R'
        dRdt = (I_n + I_s) / a_r - mu * R
        # Eq 8: D'
        dDdt = mu_h * H + mu_s * max(0, I_s - H) + mu_n * I_n
        # Eq 9: B'
        dBdt = r_b * B
        # Eq 10: C_h'
        dChdt = c_h * H
        # Eq 11: C_l'
        dCldt = c_l * (dDdt + I_n + I_s)

        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Nominal parameter values
    params = {
        'a_i': 6, 'a_r': 21, 'f_s': 0.2, 'mu': 1/45625, 'mu_n': 1/2100,
        'mu_s': 1/600, 'mu_h': 1/2400, 'beta_h': 0.1, 'c_h': 1, 'c_l': 1, 'r_b': 0
    }

    # Initial conditions
    y0 = [999999, 0, 1, 0, 0, 0, 0, 2900000, 0, 0]
    
    # Step 2: Define events and solve the ODE system
    # ------------------------------------------------
    
    # Event 1: t_n1 is when E > I_n, so we look for where E - I_n = 0
    def event_E_gt_In(t, y, params): return y[1] - y[2]
    event_E_gt_In.terminal = False
    event_E_gt_In.direction = 1 # Trigger when crossing from negative to positive

    # Event 2: t_n2 is when C_h > D, so we look for where C_h - D = 0
    def event_Ch_gt_D(t, y, params): return y[8] - y[6]
    event_Ch_gt_D.terminal = False
    event_Ch_gt_D.direction = 1

    # Time span for the simulation (in days)
    t_span = [0, 365]

    # Solve the ODE system
    sol = solve_ivp(
        fun=model,
        t_span=t_span,
        y0=y0,
        args=(params,),
        events=(event_E_gt_In, event_Ch_gt_D),
        method='RK45',
        dense_output=True
    )
    
    # Step 3: Extract results and compute final answer
    # ------------------------------------------------
    
    # Extract event times in days and convert to hours
    t_n1_days = sol.t_events[0][0]
    t_n2_days = sol.t_events[1][0]
    
    t_n1 = t_n1_days * 24  # Smallest time in hours where E > I_n
    t_n2 = t_n2_days * 24  # Smallest time in hours where C_h > D
    
    # Value of X_0 determined from the puzzle part
    X_0 = 267

    # Calculate the final answer
    final_answer = t_n2 * (X_0 - t_n1)

    # Print the components of the final calculation
    print(f"X_0 = {X_0}")
    print(f"t_n1 (hours) = {t_n1}")
    print(f"t_n2 (hours) = {t_n2}")
    print(f"Final equation: {t_n2:.4f} * ({X_0} - {t_n1:.4f})")
    print(f"Final Answer = {final_answer}")
    
solve_epidemiological_enigma()