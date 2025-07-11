import numpy as np
from scipy.integrate import solve_ivp
import math

def main():
    """
    This script solves the multi-part epidemiological modeling problem.
    """
    
    # Step 1: Define the model and parameters for simulation.
    
    def seir_model(t, y, p):
        """
        Defines the system of differential equations for the disease model.
        """
        S, E, In, Is, H, R, D, B, Ch, Cl = y

        # Quarantine-dependent contact rates
        if p["quar_start"] <= t <= p["quar_end"]:
            beta_n_t = p["beta_n_quar"]
            beta_s_t = p["beta_s_quar"]
        else:
            beta_n_t = p["beta_n_norm"]
            beta_s_t = p["beta_s_norm"]
            
        T = max(1e-9, E + In + Is + R + S) # Avoid division by zero

        # Equations as defined in the problem description
        unhosp_severe = max(0, Is - H)
        infection_rate = (beta_s_t * S * unhosp_severe / T) + \
                         (p["beta_h"] * H * S / T) + \
                         (beta_n_t * S * In / T)

        dSdt = -infection_rate - p["mu"] * S
        dEdt = infection_rate - (E * (1/p["a_i"] + p["mu"]))
        dIndt = E * (1-p["f_s"]) / p["a_i"] - In / p["a_r"] - p["mu_n"] * In
        dIsdt = E * p["f_s"] / p["a_i"] - Is / p["a_r"] - p["mu_h"] * H - p["mu_s"] * max(0, Is - H)
        
        if H < B:
            new_hosp = min(B - H, E * p["f_s"] / p["a_i"])
        else:
            new_hosp = 0
        dHdt = new_hosp - H / p["a_r"] - p["mu_h"] * H
        
        dRdt = (In + Is) / p["a_r"] - p["mu"] * R
        dDdt = p["mu_h"] * H + p["mu_s"] * max(0, Is - H) + p["mu_n"] * In
        dBdt = p["r_b"] * B
        dChdt = p["c_h"] * H
        dCldt = p["c_l"] * (dDdt + In + Is)

        return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

    nominal_params = {
        "beta_h": 0.1, "a_i": 6.0, "a_r": 21.0, "f_s": 0.2, "mu": 1.0/45625,
        "mu_n": 1.0/2100, "mu_s": 1.0/600, "mu_h": 1.0/2400, "c_h": 1.0,
        "c_l": 1.0, "r_b": 0.0, "beta_n_norm": 0.5, "beta_s_norm": 0.5,
        "beta_n_quar": 0.125, "beta_s_quar": 0.125, "quar_start": 60.0, "quar_end": 116.0
    }
    y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]

    # Simulation for high accuracy
    t_span = [0, 365]
    t_eval = np.linspace(t_span[0], t_span[1], (t_span[1] - t_span[0]) * 24 + 1)
    sol = solve_ivp(lambda t, y: seir_model(t, y, nominal_params), t_span, y0, dense_output=True, t_eval=t_eval)

    # Step 2: Find t_n1 and t_n2
    t_hours = sol.t * 24
    E_vals, In_vals, D_vals, Ch_vals = sol.y[1], sol.y[2], sol.y[6], sol.y[8]
    
    t_n1 = t_hours[np.where(E_vals > In_vals)[0][0]]
    t_n2 = t_hours[np.where(Ch_vals > D_vals)[0][0]]

    # Step 3: Use the results of the plot analysis to define p_n
    # The analysis concluded the following mapping from plot number n to parameter ID p_n
    p = {
        1: 2,   # Plot 1: S(t) vs mu_s (mortality of severe)
        2: 6,   # Plot 2: S(t) vs f_s (fraction severe)
        3: 1,   # Plot 3: S(t) vs mu (general mortality)
        4: 3,   # Plot 4: S(t) vs mu_n (mortality of normal)
        5: 7,   # Plot 5: C_l(t) vs c_l (cost of lost productivity)
        6: 5,   # Plot 6: S(t) vs a_i (incubation period)
        7: 12,  # Plot 7: C_h(t) vs c_h (healthcare cost)
        8: 8,   # Plot 8: S(t) vs mu_h (mortality of hospitalized)
        9: 9    # Plot 9: R(t) or D(t) vs beta_h (contact rate of hospitalized)
    }

    # Step 4: Calculate X_0
    X_0 = sum(n * p[n] for n in range(1, 10))

    # Step 5: Calculate the final answer and print the process
    final_answer = t_n2 * (X_0 - t_n1)

    print("--- Calculation Steps ---")
    print(f"1. Simulated model to find threshold times.")
    print(f"   Smallest time t where E > I_n (t_n1) = {t_n1:.0f} hours.")
    print(f"   Smallest time t where C_h > D (t_n2) = {t_n2:.0f} hours.")
    
    print("\n2. Deduced parameter mapping from plot analysis:")
    p_str = ', '.join([f'p_{n}={p[n]}' for n in sorted(p.keys())])
    print(f"   {p_str}")

    print("\n3. Calculated X_0 = sum(n * p_n):")
    X_0_calc_str = " + ".join([f"{n}*{p[n]}" for n in range(1, 10)])
    print(f"   X_0 = {X_0_calc_str}")
    X_0_val_str = " + ".join([str(n * p[n]) for n in range(1, 10)])
    print(f"   X_0 = {X_0_val_str} = {X_0}")

    print("\n4. Calculated final answer = t_n2 * (X_0 - t_n1):")
    print(f"   = {t_n2:.0f} * ({X_0} - {t_n1:.0f})")
    print(f"   = {t_n2:.0f} * ({X_0 - t_n1:.0f})")
    print(f"   = {final_answer:.0f}")

    print("\n--- Final Answer ---")
    print(f"{final_answer:.0f}")

if __name__ == "__main__":
    main()