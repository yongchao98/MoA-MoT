import numpy as np
from scipy.integrate import solve_ivp

def solve_epidemic_model():
    """
    This function solves the SEIR-like model for the disease outbreak,
    calculates the epidemiological thresholds t_n1 and t_n2,
    determines the parameter mapping for the puzzle to find X_0,
    and computes the final answer.
    """

    # Part 1: Solve for t_n1 and t_n2
    # ---------------------------------
    
    # Define nominal model parameters
    params = {
        'a_i': 6.0,         # Incubation period (days)
        'a_r': 21.0,        # Infectious period (days)
        'f_s': 0.2,         # Fraction developing severe symptoms
        'mu': 1.0/45625.0,  # Baseline mortality rate
        'mu_n': 1.0/2100.0, # Mortality rate for normally infected
        'mu_s': 1.0/600.0,  # Mortality rate for severely infected
        'mu_h': 1.0/2400.0, # Mortality rate for hospitalized
        'beta_h': 0.1,      # Contact rate for hospitalized
        'c_h': 1.0,         # Healthcare cost unit
        'c_l': 1.0,         # Lost productivity cost unit
        'r_b': 0.0          # Hospital bed growth rate
    }

    # Initial conditions for the state vector y = [S, E, In, Is, H, R, D, B, Ch, Cl]
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

    def get_betas(t):
        # Time-dependent contact rates
        if 60 <= t <= 116:
            return 0.125, 0.125  # Quarantine period
        return 0.5, 0.5          # Normal period

    def seir_model(t, y, p):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        beta_n_t, beta_s_t = get_betas(t)
        
        T = max(0, S + E + In + Is + R)
        if T == 0:
            return np.zeros_like(y)
            
        unhospitalized_severe = max(0, Is - H)
        force_of_infection = (beta_s_t * unhospitalized_severe + beta_n_t * In + p['beta_h'] * H) / T

        dSdt = -force_of_infection * S - p['mu'] * S
        dEdt = force_of_infection * S - E * (1/p['a_i'] + p['mu'])
        dIndt = E * (1-p['f_s'])/p['a_i'] - In/p['a_r'] - p['mu_n'] * In
        dIsdt = E * p['f_s']/p['a_i'] - Is/p['a_r'] - p['mu_h'] * H - p['mu_s'] * unhospitalized_severe
        
        if H < B:
            hosp_inflow = min(B - H, E * p['f_s'] / p['a_i'])
        else:
            hosp_inflow = 0
            
        dHdt = hosp_inflow - H/p['a_r'] - p['mu_h'] * H
        dRdt = (In + Is)/p['a_r'] - p['mu'] * R
        dDdt = p['mu_h'] * H + p['mu_s'] * unhospitalized_severe + p['mu_n'] * In
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + In + Is)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], (t_span[1] - t_span[0]) * 100 + 1)
    
    sol = solve_ivp(seir_model, t_span, y0, args=(params,), dense_output=True, t_eval=t_eval)

    # Find t_n1: smallest time t (in hours) such that E > In
    t_n1_idx = np.where(sol.y[1] > sol.y[2])[0][0]
    t_n1 = sol.t[t_n1_idx] * 24

    # Find t_n2: smallest time t (in hours) such that Ch > D
    t_n2_idx = np.where(sol.y[8] > sol.y[6])[0][0]
    t_n2 = sol.t[t_n2_idx] * 24

    # Part 2: Solve the Tracing Paper Parameter Puzzle
    # ------------------------------------------------
    # Based on qualitative analysis described in the thinking steps.
    p = {
        1: 5,   # a_i
        2: 11,  # r_b
        3: 1,   # mu
        4: 2,   # mu_s
        5: 12,  # c_h
        6: 15,  # q_f
        7: 7,   # c_l
        8: 8,   # mu_h
        9: 6    # f_s
    }

    # Part 3: Final Calculation
    # --------------------------
    X_0 = sum(n * p[n] for n in range(1, 10))
    final_answer = t_n2 * (X_0 - t_n1)

    print(f"t_n1 = {t_n1}")
    print(f"t_n2 = {t_n2}")
    print(f"X_0 = {X_0}")
    print(f"Final Answer: {t_n2} * ({X_0} - {t_n1}) = {final_answer}")
    return final_answer

final_value = solve_epidemic_model()
print(f"<<<{final_value}>>>")