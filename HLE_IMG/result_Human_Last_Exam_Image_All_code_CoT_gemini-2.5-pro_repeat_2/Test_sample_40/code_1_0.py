import math

# Define the given parameters and base values
V_initial_pu = 0.85      # Initial voltage at Bus B in per unit
V_final_pu = 1.0         # Final (target) voltage at Bus B in per unit
R_F_pu = 0.15            # Fault resistance in per unit
S_base_MVA = 100         # Base power in MVA
harmonic_loss_increase = 0.04 # 4% increase in system losses

# State the key assumption for solvability.
# The provided impedance data (Zs, Zf) is inconsistent for standard models.
# We assume an effective Thévenin reactance of the faulted system as seen from Bus B.
# A value of X_th = 0.6 p.u. is assumed as it yields a physically reasonable solution.
X_th_pu = 0.6

# --- Step 1: Calculate the Optimal Reactive Power Injection (Q_opt) ---
# We use a simplified formula relating voltage change to reactive power compensation.
# Q_opt = (V_final^2 - V_initial^2) / X_th

delta_V_squared = V_final_pu**2 - V_initial_pu**2
Q_opt_pu = delta_V_squared / X_th_pu

# Convert the reactive power from per unit to MVAR
Q_opt_MVAR = Q_opt_pu * S_base_MVA

# --- Step 2: Calculate the System's Total Real Power Losses ---
# Assume the dominant loss is the power dissipated in the fault resistance R_F.
# The voltage at the fault is assumed to be restored to 1.0 p.u. after compensation.
V_fault_pu = 1.0

# Calculate the fundamental power loss in the fault impedance
P_loss_fundamental_pu = V_fault_pu**2 / R_F_pu

# Add the 4% harmonic loss component
P_loss_total_pu = P_loss_fundamental_pu * (1 + harmonic_loss_increase)

# Convert the total power loss from per unit to MW
P_loss_total_MW = P_loss_total_pu * S_base_MVA

# --- Step 3: Print the results with detailed equations ---
print("--- Calculation of Optimal Reactive Power (Q_opt) ---")
print(f"The minimum reactive power is calculated using the equation: Q_opt = (V_final^2 - V_initial^2) / X_th")
print(f"Q_opt = ({V_final_pu}^2 - {V_initial_pu}^2) / {X_th_pu}")
print(f"Q_opt = ({V_final_pu**2} - {V_initial_pu**2:.4f}) / {X_th_pu}")
print(f"Q_opt = {delta_V_squared:.4f} / {X_th_pu}")
print(f"The optimal reactive power Q_opt is {Q_opt_pu:.4f} p.u., which is {Q_opt_MVAR:.2f} MVAR.\n")

print("--- Calculation of System Real Power Losses ---")
print("The fundamental power loss is calculated using the equation: P_loss = V_fault^2 / R_F")
print(f"P_loss_fundamental = {V_fault_pu}^2 / {R_F_pu}")
print(f"P_loss_fundamental = {V_fault_pu**2} / {R_F_pu} = {P_loss_fundamental_pu:.4f} p.u.")
print("The total power loss includes a 4% increase due to harmonics.")
print(f"P_loss_total = P_loss_fundamental * (1 + {harmonic_loss_increase})")
print(f"P_loss_total = {P_loss_fundamental_pu:.4f} * {1 + harmonic_loss_increase}")
print(f"The total system real power loss is {P_loss_total_pu:.4f} p.u., which is {P_loss_total_MW:.2f} MW.")