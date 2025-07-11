import numpy as np
import cmath

# --- 1. Given Data and Assumptions ---
# System base values
S_base = 100.0  # MVA
V_base_B = 220.0  # kV
V_base_A = 400.0  # kV

# Impedances (assuming given values are in per-unit)
Z_S_complex = complex(0.02, 0.10) # p.u.
R_S = Z_S_complex.real
X_S = Z_S_complex.imag
Z_F = 0.15 # p.u. (purely resistive)

# Voltage conditions
V_A = 1.0 # p.u. (magnitude)
V_B_pre = 0.85 # p.u. (magnitude before compensation)
V_B_post = 1.0 # p.u. (magnitude after compensation)

# STATCOM and Loss parameters
Q_max = 50.0 # MVAR
Q_max_pu = Q_max / S_base # p.u.
harmonic_loss_factor = 1.04

# --- 2. Formulate and Solve Optimization Problem ---

# The objective is to minimize Q_comp. This is achieved by maximizing the
# reactive power Q_B' delivered from the grid.
# The post-compensation power flow equation for reactive power at Bus B is:
# Q_B' = (V_A*V_B_post*(X_S*cos(d_post) - R_S*sin(d_post)) - V_B_post^2*X_S) / |Z_S|^2
# To maximize Q_B', we need to maximize the term (X_S*cos(d_post) - R_S*sin(d_post)).
# The maximum occurs when its derivative with respect to d_post is zero.
# d/d(d_post) [ ... ] = -X_S*sin(d_post) - R_S*cos(d_post) = 0
# This gives tan(d_post) = -R_S / X_S.

delta_post_rad = np.arctan(-R_S / X_S)
delta_post_deg = np.rad2deg(delta_post_rad)

# Now we have the optimal angle for the post-compensation state.
# We can use the power flow equations to solve for the unknowns.
# Let's define the power flow function for convenience.
def get_power_at_bus_b(Vb_mag, delta_b_rad):
    Vb = cmath.rect(Vb_mag, delta_b_rad)
    Z_S_conj = Z_S_complex.conjugate()
    S_b = (Vb * V_A.conjugate() - abs(Vb)**2) / Z_S_conj
    return S_b

# Post-compensation state analysis
S_B_post_grid = get_power_at_bus_b(V_B_post, delta_post_rad)
P_B_post_grid = S_B_post_grid.real
Q_B_post_grid = S_B_post_grid.imag

# The total power at Bus B is S_load + S_fault - S_comp
# S_B_post_grid = S_load + S_fault_post - j*Q_opt
P_fault_post = V_B_post**2 / Z_F
S_load_plus_P_fault_post = P_B_post_grid
P_load = S_load_plus_P_fault_post - P_fault_post

Q_load_minus_Q_opt = Q_B_post_grid

# Pre-compensation state analysis
# We need to find delta_pre, which satisfies the pre-compensation conditions
# for the determined S_load = P_load + j*Q_load.
# We have two equations for two unknowns (delta_pre, Q_load):
# P_B_pre_grid = P_load + P_fault_pre
# Q_B_pre_grid = Q_load
from scipy.optimize import fsolve

def pre_comp_equations(vars):
    q_load, delta_pre_rad = vars
    P_fault_pre = V_B_pre**2 / Z_F
    
    # Expected power at Bus B based on load and fault
    P_B_pre_expected = P_load + P_fault_pre
    Q_B_pre_expected = q_load

    # Calculated power at Bus B from power flow equation
    S_B_pre_calc = get_power_at_bus_b(V_B_pre, delta_pre_rad)
    
    eq1 = P_B_pre_expected - S_B_pre_calc.real
    eq2 = Q_B_pre_expected - S_B_pre_calc.imag
    return [eq1, eq2]

# Initial guess for the solver
initial_guess = [0.0, 0.0] 
q_load, delta_pre_rad = fsolve(pre_comp_equations, initial_guess)

# Now we have Q_load, so we can find Q_opt
# Q_load_minus_Q_opt = Q_B_post_grid => Q_opt = Q_load - Q_B_post_grid
Q_opt_pu = q_load - Q_B_post_grid
Q_opt_mvar = Q_opt_pu * S_base

# --- 3. Calculate System Losses ---
# Losses are calculated in the post-compensation state.
V_A_complex = cmath.rect(V_A, 0)
V_B_post_complex = cmath.rect(V_B_post, delta_post_rad)
I_AB_post_complex = (V_A_complex - V_B_post_complex) / Z_S_complex

P_loss_fundamental_pu = abs(I_AB_post_complex)**2 * R_S
P_loss_total_pu = P_loss_fundamental_pu * harmonic_loss_factor
P_loss_total_mw = P_loss_total_pu * S_base


# --- 4. Print Results ---
print("--- Optimization Results ---")
print(f"The minimum reactive power injection required from the STATCOM is:")
print(f"Q_opt = {Q_opt_mvar:.2f} MVAR")
print("\n--- Feasibility Check ---")
print(f"STATCOM maximum capacity: Q_max = {Q_max:.2f} MVAR")
if Q_opt_mvar > Q_max:
    print("Result is INFEASIBLE. The required reactive power exceeds the STATCOM's maximum capacity.")
else:
    print("Result is FEASIBLE. The required reactive power is within the STATCOM's capacity.")

print("\n--- System Losses Calculation ---")
print(f"The total real power losses in the system after compensation are:")
print(f"P_loss = {P_loss_total_mw:.2f} MW")

print("\n--- Detailed Calculation Steps ---")
print("1. Assumed per-unit values:")
print(f"   Z_S = {Z_S_complex.real:.2f} + j{Z_S_complex.imag:.2f} p.u.")
print(f"   Z_F = {Z_F:.2f} p.u.")
print(f"   Voltage before compensation: V_B_pre = {V_B_pre:.2f} p.u.")
print("2. Optimization to find optimal angle for voltage restoration:")
print(f"   Optimal post-compensation angle: delta_post = {delta_post_deg:.2f} degrees")
print("3. Solving for required compensation:")
print(f"   Calculated system load: P_load = {P_load:.2f} p.u., Q_load = {q_load:.2f} p.u.")
print(f"   Required compensation: Q_opt = {q_load:.2f} - ({Q_B_post_grid:.2f}) = {Q_opt_pu:.2f} p.u.")
print(f"   Equation: {Q_opt_mvar:.2f} MVAR = {Q_opt_pu:.2f} p.u. * {S_base:.1f} MVA")
print("4. Solving for system losses:")
print(f"   Post-compensation current: I_AB = {abs(I_AB_post_complex):.2f} p.u.")
print(f"   Fundamental Loss = |{abs(I_AB_post_complex):.2f}|^2 * {R_S:.2f} = {P_loss_fundamental_pu:.4f} p.u.")
print(f"   Total Loss = {P_loss_fundamental_pu:.4f} p.u. * {harmonic_loss_factor:.2f} = {P_loss_total_pu:.4f} p.u.")
print(f"   Equation: {P_loss_total_mw:.2f} MW = {P_loss_total_pu:.4f} p.u. * {S_base:.1f} MVA")
<<<182.91, 7.76>>>