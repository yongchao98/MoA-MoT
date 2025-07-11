import numpy as np
from scipy.integrate import solve_ivp

# Define model parameters
params = {
    'a_i': 6.0, 'a_r': 21.0, 'f_s': 0.2, 'mu': 1.0 / 45625.0,
    'mu_n': 1.0 / 2100.0, 'mu_s': 1.0 / 600.0, 'mu_h': 1.0 / 2400.0,
    'beta_h': 0.1, 'c_h': 1.0, 'c_l': 1.0, 'r_b': 0.0,
}

# Define time-dependent contact rates
def get_betas(t):
    return (0.125, 0.125) if 60 <= t <= 116 else (0.5, 0.5)

# Define the system of ODEs
def seir_model(t, y, p):
    S, E, I_n, I_s, H, R, D, B, C_h, C_l = y
    beta_n, beta_s = get_betas(t)
    T = max(0, E + I_n + I_s + R + S)
    if T == 0: T = 1e-9
    I_s_not_hosp = max(0, I_s - H)
    total_infection = (beta_s * S * I_s_not_hosp + p['beta_h'] * H * S + beta_n * S * I_n) / T
    dSdt = -total_infection - p['mu'] * S
    dEdt = total_infection - E * (1/p['a_i'] + p['mu'])
    dIndt = E * (1 - p['f_s']) / p['a_i'] - I_n / p['a_r'] - p['mu_n'] * I_n
    new_severe_rate = E * p['f_s'] / p['a_i']
    hospitalization_rate = min(B - H, new_severe_rate) if H < B else 0
    dHdt = hospitalization_rate - H / p['a_r'] - p['mu_h'] * H
    dIsdt = new_severe_rate - I_s / p['a_r'] - p['mu_h'] * H - p['mu_s'] * max(0, I_s - H)
    dDdt = p['mu_h'] * H + p['mu_s'] * I_s_not_hosp + p['mu_n'] * I_n
    dRdt = (I_n + I_s) / p['a_r'] - p['mu'] * R
    dBdt = p['r_b'] * B
    dChdt = p['c_h'] * H
    dCldt = p['c_l'] * (dDdt + I_n + I_s)
    return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

# Initial conditions and time span
y0 = [999999, 0, 1, 0, 0, 0, 0, 2900000, 0, 0]
t_span = [0, 365]

# Event functions to find t_n1 and t_n2
event_tn1 = lambda t, y, p: y[1] - y[2]  # E - I_n
event_tn1.direction = 1
event_tn2 = lambda t, y, p: y[8] - y[6]  # C_h - D
event_tn2.direction = 1

# Solve the ODE system
sol = solve_ivp(seir_model, t_span, y0, args=(params,), dense_output=True, events=(event_tn1, event_tn2), method='Radau')

# Extract event times in hours
tn1_hours = sol.t_events[0][0] * 24
tn2_hours = sol.t_events[1][0] * 24