import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_epidemiological_enigma():
    """
    This function solves the three-part problem by:
    1. Simulating the given ODE model to find threshold times t_n1 and t_n2.
    2. Using a pre-determined mapping of parameters to plots to calculate X_0.
    3. Calculating the final answer as specified.
    """

    # --- Part 1: ODE Simulation ---

    # Model Parameters (Nominal)
    a_i = 6.0
    a_r = 21.0
    f_s = 0.2
    mu = 1.0 / 45625.0
    mu_n = 1.0 / 2100.0
    mu_s = 1.0 / 600.0
    mu_h = 1.0 / 2400.0
    beta_h = 0.1
    c_h = 1.0
    c_l = 1.0
    r_b = 0.0

    # Quarantine parameters
    q_start = 60.0
    q_end = 116.0
    beta_normal = 0.5
    beta_quarantine = 0.125

    # Initial Conditions
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]

    # ODE System Definition
    def dydt(t, y):
        S, E, I_n, I_s, H, R, D, B, C_h, C_l = y
        
        T = max(0, S + E + I_n + I_s + R)
        if T == 0: T = 1

        beta_n = beta_quarantine if q_start <= t <= q_end else beta_normal
        beta_s = beta_n

        I_s_unhosp = max(0, I_s - H)
        infection_force = (beta_s * I_s_unhosp + beta_h * H + beta_n * I_n) * S / T

        dS_dt = -infection_force - mu * S
        dE_dt = infection_force - E * (1/a_i + mu)
        dIn_dt = E * (1 - f_s) / a_i - I_n / a_r - mu_n * I_n
        
        death_Is_H = mu_h * H + mu_s * max(0, I_s - H)
        dIs_dt = E * f_s / a_i - I_s / a_r - death_Is_H

        hosp_in_flow = min(max(0, B - H), E * f_s / a_i) if H < B else 0
        dH_dt = hosp_in_flow - H / a_r - mu_h * H
        
        dR_dt = (I_n + I_s) / a_r - mu * R
        
        dD_dt = mu_h * H + mu_s * max(0, I_s - H) + mu_n * I_n
        
        dB_dt = r_b * B
        dCh_dt = c_h * H
        dCl_dt = c_l * (dD_dt + I_n + I_s)

        return [dS_dt, dE_dt, dIn_dt, dIs_dt, dH_dt, dR_dt, dD_dt, dB_dt, dCh_dt, dCl_dt]

    # Run Simulation
    t_span = [0, 365]
    sol = solve_ivp(dydt, t_span, y0, dense_output=True, method='RK45', rtol=1e-6, atol=1e-6)

    # Find Threshold Times
    t_fine = np.linspace(t_span[0], t_span[1], 365 * 24 * 10) # 10 points per hour
    y_fine = sol.sol(t_fine)
    E_f, In_f, D_f, Ch_f = y_fine[1], y_fine[2], y_fine[6], y_fine[8]

    t_n1_raw = -1
    for i in range(len(t_fine)):
        if E_f[i] > In_f[i]:
            t_n1_raw = t_fine[i] * 24 # convert to hours
            break
    
    t_n2_raw = -1
    for i in range(1, len(t_fine)):
        if Ch_f[i] > D_f[i]:
            t_n2_raw = t_fine[i] * 24 # convert to hours
            break

    # Using ceiling as hinted by the problem's first task formulation
    t_n1 = math.ceil(t_n1_raw)
    t_n2 = math.ceil(t_n2_raw)

    # --- Part 2: Parameter Identification ---
    # The logic for this part is explained above.
    # The resulting parameter identifiers {p_n} are:
    p = {
        1: 2,   # Plot 1: S vs mu_s
        2: 9,   # Plot 2: S vs beta_h
        3: 3,   # Plot 3: S vs mu_n
        4: 15,  # Plot 4: S vs q_f
        5: 5,   # Plot 5: R vs a_i
        6: 8,   # Plot 6: S vs mu_h
        7: 7,   # Plot 7: C_l vs c_l
        8: 1,   # Plot 8: S vs mu
        9: 6    # Plot 9: I_n vs f_s
    }

    # --- Part 3: Final Calculation ---
    X_0 = sum(n * p[n] for n in range(1, 10))
    
    final_answer = t_n2 * (X_0 - t_n1)

    print(f"Threshold time t_n1 (E > I_n): {t_n1} hours")
    print(f"Threshold time t_n2 (C_h > D): {t_n2} hours")
    print("\nParameter identification results {p_n}:")
    for i in range(1,10):
        print(f"p_{i} = {p[i]}")
    
    print(f"\nCalculated X_0 = sum(n * p_n) = {X_0}")
    
    print(f"\nFinal calculation: t_n2 * (X_0 - t_n1)")
    print(f"= {t_n2} * ({X_0} - {t_n1})")
    print(f"= {t_n2} * ({X_0 - t_n1})")
    print(f"= {final_answer}")
    print(f"\n<<< {final_answer} >>>")

solve_epidemiological_enigma()