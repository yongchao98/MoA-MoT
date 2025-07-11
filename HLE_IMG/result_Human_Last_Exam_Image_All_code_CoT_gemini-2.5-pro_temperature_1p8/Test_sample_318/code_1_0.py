import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_and_print_final_equation():
    """
    This function solves the problem by:
    1. Defining the parameter-to-plot mapping and calculating X0.
    2. Setting up and solving the system of ODEs using SciPy's solve_ivp.
    3. Using the event detection feature to find t_n1 and t_n2 precisely.
    4. Printing the final equation with the calculated numerical values.
    """

    # Part 1: Tracing Paper Puzzle - Parameter Identification and X0 calculation
    # Mapping based on qualitative analysis of the ODE system's behavior.
    p_map = {
        1: 11,  # p1 -> r_b: Effect on S(t) is delayed, matches plot 1.
        2: 2,   # p2 -> mu_s: Affects T(t) with curve crossing, matches plot 2.
        3: 3,   # p3 -> mu_n: Affects S(t), standard mortality effect, matches plot 3.
        4: 6,   # p4 -> f_s: Strong effect on S(t) dynamics, matches plot 4.
        5: 12,  # p5 -> c_h: Scales C_h(t) which grows with H(t), matches plot 5.
        6: 9,   # p6 -> beta_h: Transmission from hospital, standard effect on S(t), matches plot 6.
        7: 7,   # p7 -> c_l: Scales C_l(t), which grows early, matches plot 7.
        8: 8,   # p8 -> mu_h: Mortality in hospitals, smaller effect on S(t), matches plot 8.
        9: 5,   # p9 -> a_i: Incubation period causes time-shift, matches plot 9.
    }
    X0 = sum(n * p_map[n] for n in range(1, 10))

    # Part 2: ODE solving for t_n1 and t_n2
    
    # Model Parameters
    a_i = 6.0; a_r = 21.0; f_s = 0.2; mu = 1.0/45625.0; mu_n = 1.0/2100.0
    mu_s = 1.0/600.0; mu_h = 1.0/2400.0; beta_h = 0.1; r_b = 0.0; c_h = 1.0; c_l = 1.0

    # Initial Conditions
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]

    def seir_model(t, y):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        beta_n = beta_s = 0.125 if 60.0 <= t <= 116.0 else 0.5
        
        T_pop = S + E + In + Is + R
        T = max(1e-9, T_pop)
        
        # Ensure Is >= H for stability, as H is a sub-compartment of Is
        Is_eff = max(H, Is)
        Is_minus_H = max(0, Is_eff - H)

        interaction = (beta_s*S*Is_minus_H + beta_h*H*S + beta_n*S*In) / T
        
        dSdt = -interaction - mu * S
        dEdt = interaction - E * (1/a_i + mu)
        
        flow_E_out = E / a_i
        
        dIndt = flow_E_out * (1 - f_s) - In / a_r - mu_n * In
        
        if H >= B:
            h_intake = 0
        else:
            h_intake = min(max(0, B - H), flow_E_out * f_s)
            
        dIsdt = flow_E_out*f_s - Is_eff/a_r - mu_h*H - mu_s*Is_minus_H
        dHdt = h_intake - H/a_r - mu_h*H
        dRdt = (In + Is_eff)/a_r - mu*R
        
        dDdt = mu_h*H + mu_s*Is_minus_H + mu_n*In
        dBdt = r_b * B
        dChdt = c_h * H
        dCldt = c_l * (dDdt + In + Is_eff)
        
        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Events to find t_n1 (E > In) and t_n2 (Ch > D)
    event_E_gt_In = lambda t, y: y[1] - y[2]  # E - In
    event_E_gt_In.direction = 1
    event_Ch_gt_D = lambda t, y: y[8] - y[6]  # Ch - D
    event_Ch_gt_D.direction = 1

    sol = solve_ivp(
        seir_model, [0, 365], y0, method='RK45', dense_output=True,
        events=(event_E_gt_In, event_Ch_gt_D), rtol=1e-9, atol=1e-9
    )

    t_n1_hours = sol.t_events[0][0] * 24.0 if sol.t_events[0].size > 0 else -1
    t_n2_hours = sol.t_events[1][0] * 24.0 if sol.t_events[1].size > 0 else -1

    # Part 3: Final Answer Formulation
    print(f"{t_n2_hours}*({X0}-{t_n1_hours})")

solve_and_print_final_equation()