import numpy as np
from scipy.optimize import fsolve

# --- 1. System Parameters (in per-unit) ---
# Base values
S_base = 100.0  # MVA
V_base_B = 220.0  # kV

# System impedances
# Assuming given values are in p.u. on the 100 MVA base
Z_S_complex = 0.02 + 1j * 0.10  # p.u.
R_S = Z_S_complex.real
X_S = Z_S_complex.imag
Z_S_mag_sq = np.abs(Z_S_complex)**2

Z_F = 0.15  # p.u. (Fault impedance, assumed resistive)

# Voltage conditions
V_th = 1.0  # p.u. (Thevenin voltage of the grid)
V_B_fault = 0.85  # p.u. (Voltage at Bus B during fault, before compensation)
V_B_target = 1.0  # p.u. (Target voltage at Bus B after compensation)

# Other parameters
Q_max_mvar = 50.0 # MVAR
Q_max_pu = Q_max_mvar / S_base # p.u.
harmonic_loss_increase = 0.04 # 4% increase

# --- Power flow functions ---
def p_flow(V_B, delta):
    """Calculates real power flowing into Bus B."""
    term1 = (R_S / Z_S_mag_sq) * (V_th**2 - V_th * V_B * np.cos(delta))
    term2 = (X_S / Z_S_mag_sq) * (V_th * V_B * np.sin(delta))
    return term1 + term2

def q_flow(V_B, delta):
    """Calculates reactive power flowing into Bus B."""
    term1 = (X_S / Z_S_mag_sq) * (V_th**2 - V_th * V_B * np.cos(delta))
    term2 = (R_S / Z_S_mag_sq) * (V_th * V_B * np.sin(delta))
    return term1 - term2

# --- 2. Determine Background Reactive Load (Q_L) ---
# At the fault instant (before compensation), the load at Bus B is the fault load + background load.
# Assume background real power load P_L = 0.
P_fault_load_1 = V_B_fault**2 / Z_F

# We need to find the angle delta_1 that satisfies the real power balance.
# P_flow(V_B_fault, delta_1) = P_fault_load_1
def find_delta1(delta):
    return p_flow(V_B_fault, delta) - P_fault_load_1

# Solve for delta_1
delta1_initial_guess = np.deg2rad(30) # Initial guess for the angle
delta1_solution = fsolve(find_delta1, delta1_initial_guess)[0]

# Now calculate the reactive power consumed at the bus, which is the background load Q_L.
Q_L = q_flow(V_B_fault, delta1_solution)

# --- 3. Calculate Required Reactive Power Injection (Q_opt) ---
# With compensation, V_B is restored to V_B_target (1.0 p.u.).
# The real power load is now due to the fault at the new voltage.
P_fault_load_2 = V_B_target**2 / Z_F

# Find the new angle delta_2 that satisfies the real power balance.
def find_delta2(delta):
    return p_flow(V_B_target, delta) - P_fault_load_2

# Solve for delta_2
delta2_initial_guess = np.deg2rad(40) # Initial guess for the new angle
delta2_solution = fsolve(find_delta2, delta2_initial_guess)[0]

# Calculate the new total reactive power demand at the bus.
Q_load_2 = q_flow(V_B_target, delta2_solution)

# The required STATCOM injection Q_opt is the difference between the
# background demand (Q_L) and the new demand (Q_load_2).
# Q_load_2 = Q_L - Q_opt
Q_opt_pu = Q_L - Q_load_2
Q_opt_mvar = Q_opt_pu * S_base

# --- 4. Calculate System Losses ---
# Loss = R_S * |I|^2
# I = (V_th - V_B) / Z_S
V_B_compensated = V_B_target * (np.cos(delta2_solution) - 1j * np.sin(delta2_solution))
# Note: V_B angle is negative because delta = delta_th - delta_B, and we assume delta_th = 0.
I_line = (V_th - V_B_compensated) / Z_S_complex
I_line_mag_sq = np.abs(I_line)**2

# Base losses in the line
P_loss_base_pu = R_S * I_line_mag_sq

# Total losses including harmonic effects
P_loss_total_pu = P_loss_base_pu * (1 + harmonic_loss_increase)
P_loss_total_mw = P_loss_total_pu * S_base

# --- 5. Print Results ---
print("--- Calculation Results ---")
print(f"Background reactive load (Q_L): {Q_L * S_base:.2f} MVAR")
print("\nOptimization Goal: Restore Bus B voltage to 1.0 p.u.")
print(f"Required STATCOM reactive power injection (Q_opt): {Q_opt_mvar:.2f} MVAR")
if Q_opt_mvar > Q_max_mvar:
    print(f"Note: This required injection ({Q_opt_mvar:.2f} MVAR) exceeds the STATCOM's maximum capacity of {Q_max_mvar:.2f} MVAR.")
    print("      Therefore, restoring the voltage to the nominal value is not feasible with the given STATCOM.")
else:
    print(f"This is within the STATCOM's maximum capacity of {Q_max_mvar:.2f} MVAR.")

print("\nSystem Real Power Losses (with STATCOM active):")
print(f"Total system losses: {P_loss_total_mw:.2f} MW")

# Final Answer Formatting
print("\n--- Final Answer Equations ---")
print(f"Optimal reactive power (Q_opt) = Q_L - Q_load_compensated")
print(f"Q_opt = {Q_L:.4f} p.u. - {Q_load_2:.4f} p.u. = {Q_opt_pu:.4f} p.u.")
print(f"Q_opt = {Q_opt_pu:.4f} * {S_base:.1f} MVA = {Q_opt_mvar:.2f} MVAR")
print("\n")
print(f"System Power Loss = P_loss_base * (1 + Harmonic_Increase)")
print(f"System Power Loss = {P_loss_base_pu:.4f} p.u. * (1 + {harmonic_loss_increase:.2f}) = {P_loss_total_pu:.4f} p.u.")
print(f"System Power Loss = {P_loss_total_pu:.4f} * {S_base:.1f} MVA = {P_loss_total_mw:.2f} MW")
print(f'<<<{Q_opt_mvar:.2f}, {P_loss_total_mw:.2f}>>>')
