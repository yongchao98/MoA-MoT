import numpy as np
from scipy.integrate import solve_ivp
import math

# Step 1: Define the model and parameters to find t_n1 and t_n2
# These are derived directly from the problem description.

def disease_model(t, y):
    """The SEIR model with additional compartments as described in the problem."""
    S, E, In, Is, H, R, D, B, Ch, Cl = y
    
    # Model parameters (nominal values)
    a_i = 6.0
    a_r = 21.0
    f_s = 0.2
    mu = 1.0 / 45625.0
    mu_n = 1.0 / 2100.0
    mu_s = 1.0 / 600.0
    mu_h = 1.0 / 2400.0
    beta_h = 0.1
    r_b = 0.0
    c_h = 1.0
    c_l = 1.0
    
    # Time-dependent contact rates for normally and severely symptomatic individuals
    if 60 <= t <= 116:
        beta_n, beta_s = 0.125, 0.125
    else:
        beta_n, beta_s = 0.5, 0.5
    
    # Total population for infection dynamics
    T = max(0.0, S + E + In + Is + R)
    if T == 0.0: T = 1.0 # Avoid division by zero
        
    unhospitalized_severe = max(0.0, Is - H)
    
    # Total new infections rate
    new_infections = (beta_n * S * In / T) + \
                     (beta_s * S * unhospitalized_severe / T) + \
                     (beta_h * S * H / T)
    
    # System of differential equations
    dSdt = -new_infections - mu * S
    dEdt = new_infections - E * (1.0/a_i + mu)
    dIndt = E * (1.0 - f_s) / a_i - In / a_r - mu_n * In
    dIsdt = E * f_s / a_i - Is / a_r - mu_h * H - mu_s * unhospitalized_severe
    new_Is_rate = E * f_s / a_i
    influx_H = min(B - H, new_Is_rate) if H < B else 0.0
    dHdt = influx_H - H / a_r - mu_h * H
    dRdt = (In + Is) / a_r - mu * R
    dDdt = mu_h * H + mu_s * unhospitalized_severe + mu_n * In
    dBdt = r_b * B
    dChdt = c_h * H
    dCldt = c_l * (dDdt + In + Is)
    
    return [dSdt, dEdt, dIndt, dIsdt, dHdt, dRdt, dDdt, dBdt, dChdt, dCldt]

# Initial conditions for all compartments
y0 = [999999.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2900000.0, 0.0, 0.0]

# Time span for the simulation (365 days)
t_span = [0, 365]

# Event functions to find the exact times for t_n1 and t_n2
def event_E_gt_In(t, y): return y[1] - y[2] # Condition for t_n1: E > In
event_E_gt_In.terminal = False
event_E_gt_In.direction = 1 # Find when value crosses from negative to positive

def event_Ch_gt_D(t, y): return y[8] - y[6] # Condition for t_n2: C_h > D
event_Ch_gt_D.terminal = False
event_Ch_gt_D.direction = 1

# Solve the ODE system using solve_ivp
solution = solve_ivp(disease_model, t_span, y0, method='RK45', 
                     events=[event_E_gt_In, event_Ch_gt_D], dense_output=True,
                     rtol=1e-8, atol=1e-8)

# Extract t_n1 and t_n2 from the event times and convert from days to hours
t_n1_days = solution.t_events[0][0]
t_n2_days = solution.t_events[1][0]
t_n1 = t_n1_days * 24.0
t_n2 = t_n2_days * 24.0

# Step 2: Determine the parameter mapping for each plot to find X_0
# This mapping is based on the qualitative analysis described in the thinking process.
# p_n = {p_1, p_2, ..., p_9} where p_n is the identifier for the parameter varied in plot n.
p_map = [9, 3, 1, 6, 12, 5, 7, 8, 2]

# Calculate X_0 = sum_{n=1 to 9} n * p_n
X_0 = sum((n + 1) * p_map[n] for n in range(9))

# Step 3: Calculate the final answer for the puzzle
result = t_n2 * (X_0 - t_n1)

# Print the final equation with the calculated values, as requested
print("Final Answer Calculation:")
print(f"t_n2 * (X_0 - t_n1)")
print(f"Numerically solving the system, we find:")
print(f"t_n1 (smallest time in hours where E > I_n) = {t_n1:.6f}")
print(f"t_n2 (smallest time in hours where C_h > D) = {t_n2:.6f}")
print(f"From qualitative analysis of the plots, we deduce the parameter mapping, giving:")
print(f"X_0 = {X_0}")
print(f"Plugging these values into the formula:")
print(f"{t_n2:.6f} * ({X_0} - {t_n1:.6f}) = {result:.6f}")
print(f"<<<{result:.1f}>>>")