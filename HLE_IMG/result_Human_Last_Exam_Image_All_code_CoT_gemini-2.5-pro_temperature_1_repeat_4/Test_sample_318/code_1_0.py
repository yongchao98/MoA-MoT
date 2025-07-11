import numpy as np
from scipy.integrate import solve_ivp

def solve_and_calculate():
    """
    This function solves the entire problem by first running the ODE simulation,
    then calculating the parameter puzzle result, and finally computing the final answer.
    """

    # PART 1: Epidemiological Threshold Time Calculation
    
    # Define model parameters
    params = {
        'ai': 6.0, 'ar': 21.0, 'fs': 0.2, 'mu': 1.0 / 45625.0,
        'mun': 1.0 / 2100.0, 'mus': 1.0 / 600.0, 'muh': 1.0 / 2400.0,
        'bh': 0.1, 'ch': 1.0, 'cl': 1.0, 'rb': 0.0
    }

    # Define the system of ODEs
    def model(t, y, p):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        # Time-dependent contact rates for quarantine period
        bn = bs = 0.125 if 60 <= t <= 116 else 0.5
            
        T = S + E + In + Is + R
        if T <= 0: return np.zeros_like(y)

        Is_unhosp = max(0, Is - H)
        infection_term = (bn * In + bs * Is_unhosp + p['bh'] * H) * S / T
        
        dS_dt = -infection_term - p['mu'] * S
        dE_dt = infection_term - (1/p['ai'] + p['mu']) * E
        
        from_E_to_I = E / p['ai']
        
        dIn_dt = (1 - p['fs']) * from_E_to_I - In / p['ar'] - p['mun'] * In
        dIs_dt = p['fs'] * from_E_to_I - Is / p['ar'] - p['mus'] * Is_unhosp - p['muh'] * H

        new_hosp = 0
        if B - H > 0:
            new_hosp = min(B - H, p['fs'] * from_E_to_I)
        
        dH_dt = new_hosp - H / p['ar'] - p['muh'] * H
        dR_dt = (In + Is) / p['ar'] - p['mu'] * R
        dD_dt = p['mun'] * In + p['mus'] * Is_unhosp + p['muh'] * H
        dB_dt = p['rb'] * B
        dCh_dt = p['ch'] * H
        dCl_dt = p['cl'] * (dD_dt + In + Is)
        
        return [dS_dt, dE_dt, dIn_dt, dIs_dt, dH_dt, dR_dt, dD_dt, dB_dt, dCh_dt, dCl_dt]

    # Initial conditions
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]

    # Time span for simulation
    t_span = [0, 365]
    t_eval = np.arange(t_span[0], t_span[1], 0.001) # Use a fine time step for accuracy

    # Solve the ODE system
    sol = solve_ivp(model, t_span, y0, args=(params,), dense_output=True, t_eval=t_eval)

    # Extract solution vectors
    t_sol = sol.t
    E_sol, In_sol, D_sol, Ch_sol = sol.y[1], sol.y[2], sol.y[6], sol.y[8]

    # Find tn1: smallest time t (in hours) such that E > In
    tn1_hours = -1
    for i, t in enumerate(t_sol):
        if E_sol[i] > In_sol[i]:
            tn1_hours = t * 24
            break

    # Find tn2: smallest time t (in hours) such that Ch > D
    tn2_hours = -1
    for i, t in enumerate(t_sol):
        if t > 0 and Ch_sol[i] > D_sol[i]:
            tn2_hours = t * 24
            break

    # PART 2: Tracing Paper Parameter Puzzle
    # Based on qualitative analysis of the plots and parameter effects:
    p = {
        1: 9,   # Plot 1: S vs beta_h
        2: 13,  # Plot 2: S vs q_s (quarantine start)
        3: 6,   # Plot 3: S vs f_s
        4: 3,   # Plot 4: S vs mu_n
        5: 11,  # Plot 5: R/D vs r_b
        6: 2,   # Plot 6: S vs mu_s
        7: 7,   # Plot 7: C_l vs c_l
        8: 1,   # Plot 8: T vs mu
        9: 5    # Plot 9: R vs a_i
    }
    
    # Calculate X0
    X0 = sum(n * p[n] for n in range(1, 10))

    # FINAL CALCULATION
    result = tn2_hours * (X0 - tn1_hours)

    # Print the required components of the final equation
    print(f"tn1 = {tn1_hours}")
    print(f"tn2 = {tn2_hours}")
    print(f"X0 = {X0}")
    print(f"Result = {tn2_hours} * ({X0} - {tn1_hours})")
    print(f"Final Answer = {result}")
    print(f"<<<{result}>>>")

solve_and_calculate()