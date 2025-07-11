import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_epidemiological_enigma():
    """
    Solves the multi-part epidemiological problem.
    1. Simulates the ODE model to find t_n1 and t_n2.
    2. Uses pre-determined parameter mappings to calculate X_0.
    3. Computes and prints the final answer.
    """

    # Part 1: Epidemiological Simulation
    # Initial Conditions
    y0 = [
        999999.0,  # S(0)
        0.0,       # E(0)
        1.0,       # I_n(0)
        0.0,       # I_su(0) (unhospitalized severe)
        0.0,       # H(0)
        0.0,       # R(0)
        0.0,       # D(0)
        0.0,       # C_h(0)
        0.0        # C_l(0)
    ]

    # Model Parameters
    params = {
        'a_i': 6.0, 'a_r': 21.0, 'f_s': 0.2, 'mu': 1.0 / 45625.0,
        'mu_n': 1.0 / 2100.0, 'mu_s': 1.0 / 600.0, 'mu_h': 1.0 / 2400.0,
        'c_h': 1.0, 'c_l': 1.0, 'B': 2900000.0, 'beta_h': 0.1,
        'beta_base': 0.5, 'q_s': 60.0, 'q_l': 57.0, 'q_f': 0.125
    }

    def beta_ns(t, p):
        if p['q_s'] <= t <= p['q_s'] + p['q_l'] -1:
            return p['q_f']
        else:
            return p['beta_base']

    def model(t, y, p):
        S, E, In, Isu, H, R, D, Ch, Cl = y
        
        Ist = Isu + H
        T = max(0.0, S + E + In + Ist + R)
        if T == 0: T = 1.0

        beta_val = beta_ns(t, p)
        
        # Derivatives
        interaction_Is = beta_val * S * Isu / T
        interaction_H = p['beta_h'] * H * S / T
        interaction_In = beta_val * S * In / T
        
        dSdt = -interaction_Is - interaction_H - interaction_In - p['mu'] * S
        dEdt = interaction_Is + interaction_H - interaction_In - E * (1/p['a_i'] + p['mu'])
        dIndt = E * (1 - p['f_s']) / p['a_i'] - In / p['a_r'] - p['mu_n'] * In
        
        flow_severe = E * p['f_s'] / p['a_i']
        dH_in = min(p['B'] - H, flow_severe) if H < p['B'] else 0.0
        
        dIsu_in = flow_severe - dH_in
        
        dIsudt = dIsu_in - Isu / p['a_r'] - p['mu_s'] * Isu
        dHdt = dH_in - H / p['a_r'] - p['mu_h'] * H
        
        dRdt = (In + Ist) / p['a_r'] - p['mu'] * R
        dDdt = p['mu_n'] * In + p['mu_s'] * Isu + p['mu_h'] * H
        dChdt = p['c_h'] * H
        dCldt = p['c_l'] * (dDdt + In + Ist)
        
        return [dSdt, dEdt, dIndt, dIsudt, dHdt, dRdt, dDdt, dChdt, dCldt]

    # Event functions to find t_n1 and t_n2
    def event_E_gt_In(t, y, p): return y[1] - y[2] # E - In = 0
    event_E_gt_In.direction = 1  # Find crossing from negative to positive
    event_E_gt_In.terminal = False

    def event_Ch_gt_D(t, y, p): return y[7] - y[6] # Ch - D = 0
    event_Ch_gt_D.direction = 1
    event_Ch_gt_D.terminal = False

    # Run simulation
    t_span = [0, 500]
    sol = solve_ivp(
        model, t_span, y0, args=(params,), method='RK45', dense_output=True,
        events=(event_E_gt_In, event_Ch_gt_D)
    )

    t_n1_days = sol.t_events[0][0] if len(sol.t_events[0]) > 0 else -1
    t_n2_days = sol.t_events[1][0] if len(sol.t_events[1]) > 0 else -1
    
    t_n1_hours = t_n1_days * 24
    t_n2_hours = t_n2_days * 24

    # Part 2: Parameter Puzzle Decryption
    # Based on qualitative analysis, the mapping from plot number to parameter ID is:
    # mu:1, mu_s:2, mu_n:3, a_i:5, f_s:6, c_l:7, mu_h:8, beta_h:9, r_b:11, c_h:12, q_s:13, q_l:14, q_f:15
    p = {
        1: 15,  # plot 1 -> q_f
        2: 13,  # plot 2 -> q_s
        3: 3,   # plot 3 -> mu_n
        4: 14,  # plot 4 -> q_l
        5: 9,   # plot 5 -> beta_h
        6: 2,   # plot 6 -> mu_s
        7: 7,   # plot 7 -> c_l
        8: 5,   # plot 8 -> a_i
        9: 6    # plot 9 -> f_s
    }
    
    X_0 = sum(n * p[n] for n in range(1, 10))

    # Part 3: Final Calculation
    final_answer = t_n2_hours * (X_0 - t_n1_hours)
    
    print(f"Epidemiological Threshold Time Calculation:")
    print(f"t_n1 (time when E > I_n) = {t_n1_hours:.4f} hours")
    print(f"t_n2 (time when C_h > D) = {t_n2_hours:.4f} hours")
    print(f"\nTracing Paper Parameter Puzzle Solution:")
    print(f"p = {list(p.values())}")
    print(f"X_0 = sum(n * p_n) = {X_0}")
    print(f"\nThe Epidemiological Enigma's Cipher:")
    print(f"Result = t_n2 * (X_0 - t_n1)")
    print(f"Result = {t_n2_hours:.4f} * ({X_0} - {t_n1_hours:.4f})")
    print(f"Result = {final_answer:.4f}")
    
    # Final answer in specified format
    print(f"\n<<<Result: {final_answer}>>>")

solve_epidemiological_enigma()