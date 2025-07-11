import cmath
import math

# --- Given Parameters ---
V_B_sag = 0.85  # Voltage at Bus B during sag (pu)
V_B_final = 1.0  # Target restored voltage at Bus B (pu)
R_th = 0.02  # Assumed Thevenin resistance (pu)
X_th = 0.10  # Assumed Thevenin reactance (pu)
S_base = 100.0  # Base power (MVA)
harmonic_loss_increase = 0.04  # 4% increase in losses

# --- Initial Calculations ---
V_th = 1.0  # Thevenin voltage of the grid (pu)
Z_th_sq = R_th**2 + X_th**2

# Step 1: Determine the purely resistive load P_L that causes the sag to V_B = 0.85 pu
# We use the voltage-power equation with V_B = 0.85, Q_L = 0
# |V_B|^4 + |V_B|^2 * (2*P_L*R_th - |V_th|^2) + P_L^2 * |Z_th|^2 = 0
# This is a quadratic equation for P_L: a*P_L^2 + b*P_L + c = 0
a1 = Z_th_sq
b1 = 2 * R_th * V_B_sag**2
c1 = V_B_sag**4 - V_B_sag**2 * V_th**2

# Solve the quadratic equation for P_L
discriminant1 = b1**2 - 4 * a1 * c1
P_L1 = (-b1 + math.sqrt(discriminant1)) / (2 * a1)
P_L2 = (-b1 - math.sqrt(discriminant1)) / (2 * a1)
# We choose the positive power solution
P_L = P_L1 if P_L1 > 0 else P_L2
Q_L = 0.0

print("--- Step 1: Determine Load Characteristics ---")
print("To minimize reactive power injection, we assume a purely resistive load (PF=1).")
print("We find the load P_L that causes voltage to drop to 0.85 pu.")
print("The governing equation is a*P_L^2 + b*P_L + c = 0, where:")
print(f"a = |Z_th|^2 = {a1:.4f}")
print(f"b = 2*R_th*V_B_sag^2 = {b1:.4f}")
print(f"c = V_B_sag^4 - V_B_sag^2 = {c1:.4f}")
print(f"Solving the equation: {a1:.4f}*P_L^2 + {b1:.4f}*P_L + {c1:.4f} = 0")
print(f"The resulting real power load is P_L = {P_L:.4f} pu\n")


# Step 2: Determine the required reactive power Q_opt
# Now we set V_B = 1.0 and solve for the new net reactive power Q_B at the bus.
# The real power load P_B is still P_L. The new reactive load is Q_B = Q_L - Q_opt = -Q_opt
# |V_B|^4 + |V_B|^2 * (2*(P_L*R_th + Q_B*X_th) - |V_th|^2) + (P_L^2 + Q_B^2)*|Z_th|^2 = 0
# This is a quadratic equation for Q_B: a*Q_B^2 + b*Q_B + c = 0
a2 = Z_th_sq * V_B_final**2
b2 = 2 * X_th * V_B_final**2
c2 = V_B_final**4 + V_B_final**2 * (2 * P_L * R_th - V_th**2) + P_L**2 * Z_th_sq

# Solve the quadratic equation for Q_B
discriminant2 = b2**2 - 4 * a2 * c2
# We need to select the correct solution. The stable operating point corresponds to a smaller reactive power load.
Q_B1 = (-b2 + math.sqrt(discriminant2)) / (2 * a2)
Q_B2 = (-b2 - math.sqrt(discriminant2)) / (2 * a2)
# Choose the solution with the smaller magnitude, which represents the stable high-voltage solution
Q_B_final = Q_B1 if abs(Q_B1) < abs(Q_B2) else Q_B2

# Q_B_final = Q_L - Q_opt. Since Q_L = 0, Q_opt = -Q_B_final
Q_opt_pu = -Q_B_final
Q_opt_MVAR = Q_opt_pu * S_base

print("--- Step 2: Calculate Required Reactive Power (Q_opt) ---")
print("With P_L found, we calculate the net reactive power Q_B required to achieve V_B = 1.0 pu.")
print("The governing equation is a*Q_B^2 + b*Q_B + c = 0, where:")
print(f"a = |Z_th|^2 = {a2:.4f}")
print(f"b = 2*X_th = {b2:.4f}")
print(f"c = 1 + 2*P_L*R_th - 1 + P_L^2*|Z_th|^2 = {c2:.4f}")
print(f"Solving the equation: {a2:.4f}*Q_B^2 + {b2:.4f}*Q_B + {c2:.4f} = 0")
print(f"The post-compensation reactive power at the bus is Q_B = {Q_B_final:.4f} pu")
print("\nThe optimal reactive power from the STATCOM is Q_opt = -Q_B:")
print(f"Q_opt = {Q_opt_pu:.4f} pu")
print(f"Q_opt = {Q_opt_MVAR:.2f} MVAR\n")

# Step 3: Calculate System Losses
# P_loss = |I|^2 * R_th = (P_B^2 + Q_B^2) / |V_B|^2 * R_th
P_loss_base_pu = (P_L**2 + Q_B_final**2) / V_B_final**2 * R_th
P_loss_final_pu = P_loss_base_pu * (1 + harmonic_loss_increase)
P_loss_final_MW = P_loss_final_pu * S_base

print("--- Step 3: Calculate System Real Power Losses ---")
print("The real power loss in the system is calculated for the final compensated state.")
print("Base Loss Equation: P_loss_base = (P_L^2 + Q_B^2) / V_B_final^2 * R_th")
print(f"P_loss_base = (({P_L:.4f})^2 + ({Q_B_final:.4f})^2) / ({V_B_final:.1f})^2 * {R_th:.2f} = {P_loss_base_pu:.4f} pu")
print(f"Accounting for the {harmonic_loss_increase*100}% increase due to harmonics:")
print("Final Loss Equation: P_loss_final = P_loss_base * (1 + Harmonic_Increase)")
print(f"P_loss_final = {P_loss_base_pu:.4f} * {1 + harmonic_loss_increase} = {P_loss_final_pu:.4f} pu")
print(f"P_loss_final = {P_loss_final_MW:.2f} MW\n")

print("--- Final Answer ---")
print(f"Optimal Reactive Power Injection (Q_opt): {Q_opt_MVAR:.2f} MVAR")
print(f"Total System Real Power Losses: {P_loss_final_MW:.2f} MW")
<<<Optimal Reactive Power Injection (Q_opt): 126.43 MVAR, Total System Real Power Losses: 24.84 MW>>>