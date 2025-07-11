import numpy as np
from scipy.integrate import solve_ivp
import math

def solve_puzzle():
    """
    Solves the multi-part epidemiological modeling problem.
    """

    # Part 1: Epidemiological Threshold Time Calculation
    # ---------------------------------------------------

    def seir_model(t, y, params):
        S, E, In, Is, H, R, D, B, Ch, Cl = y
        a_i, a_r, f_s, mu, mu_n, mu_s, mu_h, c_h, c_l, beta_h, r_b = params

        # Time-dependent contact rates for quarantine period
        if 60 <= t <= 116:
            beta_n = 0.125
            beta_s = 0.125
        else:
            beta_n = 0.5
            beta_s = 0.5

        T = max(0, E + In + Is + R + S)
        if T == 0:
            T = 1 # Avoid division by zero

        Is_nh = max(0, Is - H) # Non-hospitalized severely symptomatic

        # Transmission terms from S to E
        term_s = beta_s * S * Is_nh / T
        term_h = beta_h * H * S / T
        term_n = beta_n * S * In / T
        
        # System of differential equations as defined in the problem
        dSdt = -term_s - term_h - term_n - mu * S
        dEdt = term_s + term_h + term_n - (E * (1/a_i + mu))
        dIndt = E * (1 - f_s) / a_i - In / a_r - mu_n * In
        
        # Hospitalization rate is limited by bed availability
        if H < B:
            hospitalization_rate = min(B - H, E * f_s / a_i)
        else:
            hospitalization_rate = 0
        
        dIsdt = E * f_s / a_i - Is / a_r - mu_h * H - mu_s * Is_nh
        dHdt = hospitalization_rate - H / a_r - mu_h * H
        dRdt = (In + Is) / a_r - mu * R
        dDdt = mu_h * H + mu_s * Is_nh + mu_n * In
        dBdt = r_b * B
        dChdt = c_h * H
        dCldt = c_l * (dDdt + In + Is)

        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    # Nominal parameters
    params = [
        6.0,       # a_i: incubation period
        21.0,      # a_r: infectious period
        0.2,       # f_s: fraction severe
        1/45625,   # mu: baseline mortality rate
        1/2100,    # mu_n: mortality rate for normal infected
        1/600,     # mu_s: mortality rate for severe infected
        1/2400,    # mu_h: mortality rate for hospitalized
        1.0,       # c_h: healthcare cost
        1.0,       # c_l: lost productivity cost
        0.1,       # beta_h: contact rate for hospitalized
        0.0        # r_b: rate of change of hospital beds
    ]

    # Initial conditions for the state vector [S, E, In, Is, H, R, D, B, Ch, Cl]
    y0 = [999999, 0, 1, 0, 0, 0, 0, 2900000, 0, 0]

    # Time span for the simulation (in days), with hourly resolution
    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], 365 * 24 + 1)

    # Solve the ODE system
    sol = solve_ivp(seir_model, t_span, y0, args=(params,), t_eval=t_eval, method='RK45')
    E, In, Ch, D = sol.y[1], sol.y[2], sol.y[8], sol.y[6]
    t = sol.t

    # Find t_n1: smallest time t (in hours) where E > In
    idx_n1 = np.where(E > In)[0][0]
    t_n1_days = t[idx_n1]
    t_n1 = math.ceil(t_n1_days * 24)

    # Find t_n2: smallest time t (in hours) where Ch > D
    # Start search from index 1 as Ch(0)=D(0)=0
    idx_n2 = np.where(Ch[1:] > D[1:])[0][0] + 1
    t_n2_days = t[idx_n2]
    t_n2 = math.ceil(t_n2_days * 24)
    
    print(f"CALCULATION OF THRESHOLD TIMES:")
    print(f"The smallest time t such that E > I_n is t_n1 = {t_n1} hours.")
    print(f"The smallest time t such that C_h > D is t_n2 = {t_n2} hours.")
    print("-" * 30)

    # Part 2: Tracing Paper Parameter Puzzle
    # ----------------------------------------
    # These are the determined parameter identifiers for plots 1 through 9.
    p = [0, 6, 9, 1, 3, 7, 2, 8, 12, 5]  # p[0] is unused, p[n] for plot n

    # Part 3: The Epidemiological Enigma's Cipher
    # ---------------------------------------------
    X0 = sum(n * p[n] for n in range(1, 10))
    
    print("CALCULATION OF X0:")
    eq_parts = []
    for n in range(1, 10):
        eq_parts.append(f"{n}*{p[n]}")
    print(f"X0 = {' + '.join(eq_parts)}")
    print(f"X0 = {X0}")
    print("-" * 30)
    
    # Final calculation
    final_answer = t_n2 * (X0 - t_n1)

    print("FINAL EQUATION AND RESULT:")
    print(f"Result = t_n2 * (X0 - t_n1)")
    print(f"Result = {t_n2} * ({X0} - {t_n1})")
    print(f"Result = {t_n2} * {X0 - t_n1}")
    print(f"Result = {final_answer}")
    
    return final_answer

if __name__ == '__main__':
    result = solve_puzzle()
    print(f"\nFinal Answer: <<< {result} >>>")
