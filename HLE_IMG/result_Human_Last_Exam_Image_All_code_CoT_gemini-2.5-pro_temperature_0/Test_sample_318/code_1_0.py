import numpy as np
from scipy.integrate import solve_ivp

def solve_epidemiological_enigma():
    """
    This function solves the entire problem by first simulating the ODE model to find t_n1 and t_n2,
    then calculating X_0 based on the solved parameter puzzle, and finally computing the final answer.
    """

    # --- Part 1: Numerical Simulation ---

    # 1.1. Define Model Parameters and Initial Conditions
    params = {
        'a_i': 6.0,       # incubation period
        'a_r': 21.0,      # infectious period
        'f_s': 0.2,       # fraction severe
        'mu': 1.0/45625,  # baseline mortality
        'mu_n': 1.0/2100, # normal infected mortality
        'mu_s': 1.0/600,  # severe infected mortality
        'mu_h': 1.0/2400, # hospitalized mortality
        'beta_h': 0.1,    # contact rate hospitalized
        'c_h': 1.0,       # healthcare cost rate
        'c_l': 1.0,       # lost productivity cost rate
        'r_b': 0.0,       # rate of change of hospital beds
    }

    # Initial state vector: [S, E, In, Is, H, R, D, B, Ch, Cl]
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]

    # 1.2. Define Time-Dependent Contact Rates
    def beta_n(t):
        return 0.125 if 60 <= t <= 116 else 0.5

    def beta_s(t):
        return 0.125 if 60 <= t <= 116 else 0.5

    # 1.3. Define the System of ODEs
    def model(t, y, p):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        
        T = max(0, E + In + Is + R + S)
        if T == 0: T = 1

        Is_nh = max(0, Is - H)
        
        transmission = (beta_s(t) * S * Is_nh + p['beta_h'] * H * S + beta_n(t) * S * In) / T
        
        dSdt = -transmission - p['mu'] * S
        dEdt = transmission - E * (1/p['a_i'] + p['mu'])
        dIndt = E * (1 - p['f_s']) / p['a_i'] - In / p['a_r'] - p['mu_n'] * In
        
        # I_s is total severe cases, H is hospitalized subset. This formulation is consistent.
        dIsdt = E * p['f_s'] / p['a_i'] - Is / p['a_r'] - p['mu_h'] * H - p['mu_s'] * Is_nh
        
        new_hosp = 0
        if H < B:
            # Rate of new hospitalizations is limited by bed availability and new severe cases
            new_hosp = min(B - H, E * p['f_s'] / p['a_i'])
        dHdt = new_hosp - H / p['a_r'] - p['mu_h'] * H
        
        dRdt = (In + Is) / p['a_r'] - p['mu'] * R
        dDdt = p['mu_h'] * H + p['mu_s'] * Is_nh + p['mu_n'] * In
        dBdt = p['r_b'] * B
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + In + Is)

        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # 1.4. Solve the ODE System
    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], t_span[1] * 24 + 1) # Hourly evaluation

    sol = solve_ivp(lambda t, y: model(t, y, params), t_span, y0, t_eval=t_eval, method='RK45')

    # 1.5. Find Threshold Times t_n1 and t_n2
    E_sol, In_sol, D_sol, Ch_sol = sol.y[1], sol.y[2], sol.y[6], sol.y[8]
    
    t_n1_hours = -1
    for i in range(len(sol.t)):
        if E_sol[i] > In_sol[i]:
            t_n1_hours = sol.t[i] * 24
            break
            
    t_n2_hours = -1
    for i in range(len(sol.t)):
        if Ch_sol[i] > 0 and D_sol[i] > 0 and Ch_sol[i] > D_sol[i]:
            t_n2_hours = sol.t[i] * 24
            break

    # --- Part 2: Calculate X_0 from Puzzle Solution ---
    # Parameter-identifier mapping:
    # mu:1, mu_s:2, mu_n:3, a_i:5, f_s:6, c_l:7, mu_h:8, beta_h:9, r_b:11, c_h:12, q_s:13, q_l:14, q_f:15
    # Based on qualitative analysis, the solution to the puzzle {p1, ..., p9} is:
    p = [
        11,  # p1: Plot 1 -> r_b
        9,   # p2: Plot 2 -> beta_h
        3,   # p3: Plot 3 -> mu_n
        15,  # p4: Plot 4 -> q_f
        6,   # p5: Plot 5 -> f_s
        2,   # p6: Plot 6 -> mu_s
        7,   # p7: Plot 7 -> c_l
        8,   # p8: Plot 8 -> mu_h
        5    # p9: Plot 9 -> a_i
    ]

    X_0 = sum((n + 1) * p[n] for n in range(9))

    # --- Part 3: Final Calculation ---
    result = t_n2_hours * (X_0 - t_n1_hours)

    # Print the components of the final equation as requested
    print(f"The final calculation is based on the formula: t_n2 * (X_0 - t_n1)")
    print(f"t_n2 = {t_n2_hours}")
    print(f"X_0 = {X_0}")
    print(f"t_n1 = {t_n1_hours}")
    print(f"Final Equation: {t_n2_hours} * ({X_0} - {t_n1_hours}) = {result}")
    print(f"<<<{result}>>>")

solve_epidemiological_enigma()