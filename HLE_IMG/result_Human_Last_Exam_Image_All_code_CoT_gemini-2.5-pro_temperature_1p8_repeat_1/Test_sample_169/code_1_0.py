import math

# --- 1. System Parameters ---
# The target voltage at Bus-W
V_target = 0.575  # in kV
# Transmission impedance from Bus-W to Fault
Z_WF_real = 0.01   # in Ohms
Z_WF_imag = 0.05   # in Ohms
# Fault resistance
R_F = 0.1        # in Ohms
# Power factor constraint for the generator
PF_min = 0.95      # lagging
# Harmonic loss factor (2% + 4% = 6% increase)
k_loss = 1.06

# --- 2. Fault Analysis Calculation ---
# Total impedance of the fault path
Z_total_real = Z_WF_real + R_F
Z_total_imag = Z_WF_imag
# Squared magnitude of the total impedance
Z_total_mag_sq = Z_total_real**2 + Z_total_imag**2

# Power is calculated as S = |V|^2 / Z_conj = |V|^2 * Z / |Z|^2
# |V_target|^2
V_target_sq = V_target**2

# Calculate the fundamental active and reactive power drawn by the fault
# P = |V|^2 * R / |Z|^2 (in MW since V is in kV)
P_WF = (V_target_sq * Z_total_real) / Z_total_mag_sq
# Q = |V|^2 * X / |Z|^2 (in MVAR since V is in kV)
Q_WF = (V_target_sq * Z_total_imag) / Z_total_mag_sq

# --- 3. Generator Power Calculation with Losses ---
# The total active power that must be generated includes harmonic losses
P_gen = P_WF * k_loss

# --- 4. Apply Generator Power Factor Constraint ---
# To minimize Q_comp, we must maximize Q_gen.
# Q_gen is maximized when the generator operates at the minimum power factor of 0.95.
# Get the angle phi from the power factor, PF = cos(phi)
phi_max_rad = math.acos(PF_min)
# The maximum reactive power from the generator is P_gen * tan(phi)
Q_gen_max = P_gen * math.tan(phi_max_rad)

# --- 5. Calculate Optimal Reactive Power Injection ---
# The reactive power from the compensator makes up the difference.
# Q_WF = Q_gen + Q_comp  =>  Q_comp = Q_WF - Q_gen
Q_opt = Q_WF - Q_gen_max

# --- 6. Output the Results ---
print("--- Problem Formulation and Solution ---")
print("Objective: Minimize Q_comp")
print("Constraint 1: Restore Bus-W voltage to {:.3f} kV".format(V_target))
print("Constraint 2: System power factor >= {:.2f} lagging".format(PF_min))
print("\n--- Calculation Steps ---")
print("1. Reactive power required by the fault (Q_WF):")
print("   Q_WF = (|V_target|^2 * X_total) / |Z_total|^2")
print("   Q_WF = ({:.3f}^2 kV^2 * {:.2f} Ohm) / {:.4f} Ohm^2 = {:.4f} MVAR".format(V_target, Z_total_imag, Z_total_mag_sq, Q_WF))
print("\n2. Max reactive power from generator (Q_gen_max) at PF=0.95:")
print("   P_gen = P_WF * loss_factor = ({:.4f} MW * {:.2f}) = {:.4f} MW".format(P_WF, k_loss, P_gen))
print("   Q_gen_max = P_gen * tan(acos({:.2f}))".format(PF_min))
print("   Q_gen_max = {:.4f} MW * {:.4f} = {:.4f} MVAR".format(P_gen, math.tan(phi_max_rad), Q_gen_max))

print("\n--- Final Equation for Optimal Reactive Power ---")
print("The optimal reactive power Q_opt is the difference between the required and generated reactive power.")
# The final requested format: output each number in the final equation
print("\nFinal Equation:")
print("Q_opt = Q_WF - Q_gen_max")
print("Q_opt = {:.4f} MVAR - {:.4f} MVAR".format(Q_WF, Q_gen_max))

print("\n--- Result ---")
print("The optimal reactive power injection required from the compensating device is: {:.4f} MVAR".format(Q_opt))

# The final answer in the requested format
final_answer = round(Q_opt, 4)
# print(f"\n<<<{final_answer}>>>")