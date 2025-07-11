import cmath
import math

# Step 1 & 2: Define constants and convert to per-unit (p.u.)
# The problem statement has inconsistencies in voltage levels. We will assume the
# numerical values are given in a consistent p.u. system, with V_target=0.575 p.u.
# and S_base = 100 MVA.
S_base_MVA = 100.0

# Given values (assumed to be in p.u.)
R_WF_pu = 0.01
X_WF_pu = 0.05
Z_WF_pu = complex(R_WF_pu, X_WF_pu)
R_F_pu = 0.1
V_W_target_pu = 0.575
Q_max_MVAR = 10.0
Q_max_pu = Q_max_MVAR / S_base_MVA
PF_min = 0.95
harmonic_loss_factor = 1.06 # 1 + 2% + 4%

# Step 3: Formulate the electrical model equations
# Model: Source at V_W, connected to Z_WF in series with R_F.
# The voltage at the fault point F is V_F.
# V_F = V_W - I * Z_WF
# V_F = I * R_F
# Solving for V_F: V_W - I * Z_WF = I * R_F => V_W = I * (R_F + Z_WF)
# So, the current I is I = V_W / (R_F + Z_WF)
# And the voltage at the fault is V_F = (V_W / (R_F + Z_WF)) * R_F

# The complex power S_W supplied at Bus-W is V_W * conj(I)
# S_W = V_W * conj(V_W / (R_F + Z_WF)) = |V_W|^2 / conj(R_F + Z_WF)

# Step 4: Incorporate harmonic losses.
# The total active power P_W supplied must cover the loss in Z_WF and the power dissipated in R_F.
# The harmonic losses specifically increase the resistive losses in the line.
# P_W = P_loss_in_Z_WF_total + P_dissipated_in_RF
# P_loss_in_Z_WF_total = harmonic_loss_factor * |I|^2 * R_WF_pu
# P_dissipated_in_RF = |I|^2 * R_F_pu
# P_W = |I|^2 * (harmonic_loss_factor * R_WF_pu + R_F_pu)
# The net reactive power Q_net = |I|^2 * X_WF_pu
#
# Let's calculate I first using V_W_target_pu
Z_total_for_current = R_F_pu + Z_WF_pu
I_W_ phasor = V_W_target_pu / Z_total_for_current
I_W_mag_sq = abs(I_W_phasor)**2

# Now calculate P_W and Q_net
P_W_pu = I_W_mag_sq * (harmonic_loss_factor * R_WF_pu + R_F_pu)
Q_net_pu = I_W_mag_sq * X_WF_pu

# Step 5: Solve the optimization problem
# Objective: minimize Q_comp_pu
# Power Balance: Q_net_pu = Q_W_pu + Q_comp_pu (Generator and Compensator supply the load)
# So, Q_comp_pu = Q_net_pu - Q_W_pu
# To minimize Q_comp_pu, we must maximize Q_W_pu.

# Constraint: Power factor >= 0.95 lagging for the generator
# This means Q_W_pu / P_W_pu <= tan(acos(PF_min))
max_Q_over_P_ratio = math.tan(math.acos(PF_min))
Q_W_max_pu = P_W_pu * max_Q_over_P_ratio

# Calculate the optimal (required) reactive power injection
Q_opt_pu = Q_net_pu - Q_W_max_pu

# Convert optimal Q back to MVAR
Q_opt_MVAR = Q_opt_pu * S_base_MVA

# Display the results step-by-step
print("--- System Parameters (p.u.) ---")
print(f"Target Voltage at Bus-W (|V_W|): {V_W_target_pu:.4f} p.u.")
print(f"Line Impedance (Z_WF): {Z_WF_pu.real:.4f} + j{Z_WF_pu.imag:.4f} p.u.")
print(f"Fault Resistance (R_F): {R_F_pu:.4f} p.u.")
print(f"Power Factor Constraint (PF): >= {PF_min}")
print(f"Harmonic Loss Factor: {harmonic_loss_factor}")
print("\n--- Calculation Steps ---")
print(f"Total Current Magnitude Squared (|I_W|^2): {I_W_mag_sq:.4f} p.u.")
print(f"Total Active Power Demand at Bus-W (P_W): {P_W_pu:.4f} p.u.")
print(f"Net Reactive Power Demand at Bus-W (Q_net): {Q_net_pu:.4f} p.u.")
print("\n--- Optimization ---")
print("Objective: Minimize Q_comp, where Q_comp = Q_net - Q_W")
print("This requires maximizing the generator's reactive power output, Q_W.")
print(f"Maximum Q_W based on PF constraint (Q_W_max): {P_W_pu:.4f} * tan(acos({PF_min})) = {Q_W_max_pu:.4f} p.u.")
print("\n--- Final Equation and Solution ---")
# Check if the required compensation is within the device's limits
if Q_opt_pu < 0:
    print("The system is too reactive. The required compensation is negative (absorption), which is outside the scope of the device.")
    print(f"Required Q_comp = Q_net - Q_W_max = {Q_net_pu:.4f} - {Q_W_max_pu:.4f} = {Q_opt_pu:.4f} p.u.")
elif Q_opt_pu > Q_max_pu:
    print("Warning: Required reactive power exceeds the compensator's maximum capacity.")
    print(f"Required Q_comp ({Q_opt_pu:.4f} p.u.) > Q_max ({Q_max_pu:.4f} p.u.).")
    print("The target voltage cannot be reached under the specified constraints.")
    print("The optimal feasible injection is the maximum capacity.")
    Q_opt_pu = Q_max_pu
    Q_opt_MVAR = Q_max_MVAR
    print(f"Optimal Feasible Q_comp = {Q_opt_pu:.4f} p.u.")

print("\nFinal Equation for Optimal Reactive Power (Q_opt):")
print(f"Q_opt = Q_net - Q_W_max")
print(f"Q_opt = {Q_net_pu:.4f} - {Q_W_max_pu:.4f}")
print(f"Q_opt = {Q_opt_pu:.4f} p.u.")
print(f"Q_opt = {Q_opt_MVAR:.4f} MVAR")

# The final answer in the requested format
final_answer = Q_opt_MVAR
# print(f"\n<<<{final_answer:.4f}>>>")