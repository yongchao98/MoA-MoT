import numpy as np
from scipy.integrate import solve_ivp
import math

# Step 1: Numerical Simulation to find t_n1 and t_n2

# Model parameters for the first part of the problem
params = {
    'a_i': 6,
    'a_r': 21,
    'f_s': 0.2,
    'mu': 1 / 45625,
    'mu_n': 1 / 2100,
    'mu_s': 1 / 600,
    'mu_h': 1 / 2400,
    'beta_h': 0.1,
    'r_b': 0,
    'c_h': 1,
    'c_l': 1
}

# Initial conditions
y0 = [
    999999,  # S(0)
    0,       # E(0)
    1,       # I_n(0)
    0,       # I_s(0)
    0,       # H(0)
    0,       # R(0)
    0,       # D(0)
    2900000, # B(0)
    0,       # C_h(0)
    0        # C_l(0)
]

def beta_ns(t):
    """ Time-dependent contact rate for non-hospitalized individuals """
    if 60 <= t <= 116:
        return 0.125  # Quarantine period
    else:
        return 0.5

def model(t, y, p):
    """ System of ODEs for the epidemiological model """
    S, E, I_n, I_s, H, R, D, B, C_h, C_l = y
    
    T = max(0.1, S + E + I_n + I_s + R) # Total living population, avoid division by zero
    
    bn_t = beta_ns(t)
    bs_t = beta_ns(t)
    
    Is_minus_H = max(0, I_s - H)
    
    infection_term = (bs_t * Is_minus_H + p['beta_h'] * H + bn_t * I_n) * S / T
    
    dSdt = -infection_term - p['mu'] * S
    dEdt = infection_term - E / p['a_i'] - p['mu'] * E
    dIndt = E * (1 - p['f_s']) / p['a_i'] - I_n / p['a_r'] - p['mu_n'] * I_n
    dIsdt = E * p['f_s'] / p['a_i'] - I_s / p['a_r'] - p['mu_h'] * H - p['mu_s'] * Is_minus_H
    
    if H < B:
        new_hospitalizations = min(B - H, E * p['f_s'] / p['a_i'])
    else:
        new_hospitalizations = 0
    dHdt = new_hospitalizations - H / p['a_r'] - p['mu_h'] * H
    
    dRdt = (I_n + I_s) / p['a_r'] - p['mu'] * R
    dDdt = p['mu_h'] * H + p['mu_s'] * Is_minus_H + p['mu_n'] * I_n
    dBdt = p['r_b'] * B
    dChdt = p['c_h'] * H
    dCldt = p['c_l'] * (dDdt + I_n + I_s)
    
    return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

# Event functions to find t_n1 and t_n2
def event_e_gt_in(t, y, p): return y[1] - y[2]  # E - I_n
event_e_gt_in.direction = 1  # Trigger when going from negative to positive
def event_ch_gt_d(t, y, p): return y[8] - y[6]  # C_h - D
event_ch_gt_d.direction = 1

t_span = [0, 365]
solution = solve_ivp(model, t_span, y0, args=(params,), dense_output=True,
                     events=(event_e_gt_in, event_ch_gt_d), max_step=0.1)

# Extract event times
t_n1_days = solution.t_events[0][0]
t_n2_days = solution.t_events[1][0]

# Convert to hours
t_n1 = t_n1_days * 24
t_n2 = t_n2_days * 24

# Step 2: Parameter identification result
# Based on qualitative analysis of the plots, as described in the thinking process.
p_identifiers = [2, 6, 9, 3, 8, 13, 7, 1, 5]

# Step 3: Final Calculation
X0 = sum((n + 1) * p for n, p in enumerate(p_identifiers))
final_answer = t_n2 * (X0 - t_n1)

print("--- Epidemiological Threshold Times ---")
print(f"t_n1 (first time E > I_n): {t_n1:.4f} hours")
print(f"t_n2 (first time C_h > D): {t_n2:.4f} hours")
print("\n--- Cipher Calculation ---")
print(f"Sum of weighted parameter identifiers (X0): {X0}")
print("\n--- Final Answer Calculation ---")
print(f"Formula: t_n2 * (X0 - t_n1)")
print(f"Values: {t_n2:.4f} * ({X0} - {t_n1:.4f})")
print(f"Result: {final_answer:.4f}")
print(f"<<<{final_answer}>>>")