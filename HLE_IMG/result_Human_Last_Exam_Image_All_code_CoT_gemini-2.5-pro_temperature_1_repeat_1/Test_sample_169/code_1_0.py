import cmath
import math

# --- Step 1: Define System Parameters in Per-Unit (p.u.) ---
# Base Power
S_base_MVA = 100.0

# Voltage at Bus-W (target)
V_W_pu = 0.575

# Transmission Impedance (assuming given values are in p.u.)
R_WF_pu = 0.01
X_WF_pu = 0.05
Z_WF_pu = complex(R_WF_pu, X_WF_pu)

# Fault Resistance (assuming given value is in p.u.)
R_F_pu = 0.1

# Maximum reactive power capacity of the compensator
Q_max_MVAR = 10.0
Q_max_pu = Q_max_MVAR / S_base_MVA

# Power factor constraint
PF_min = 0.95
# tan(phi) = Q/P. For a given PF, |Q/P| <= tan(acos(PF))
tan_phi_max = math.tan(math.acos(PF_min))

# Harmonic loss factor (6% increase in active power losses)
loss_factor = 1.06

# --- Step 2 & 3: Calculate Power Drawn by the Faulty Line ---
# Simplified fault model: Z_fault = Z_WF + R_F
Z_fault_pu = Z_WF_pu + R_F_pu

# Power S_W = V_W^2 / Z_fault*
S_W_pu = (V_W_pu**2) / Z_fault_pu.conjugate()
P_W_pu = S_W_pu.real
Q_W_pu = S_W_pu.imag

# --- Step 4: Incorporate Harmonic Losses ---
# The active power required from the generator is higher due to harmonic losses
P_gen_pu = P_W_pu * loss_factor

# --- Step 5 & 6: Apply Constraints and Optimize ---
# From PF constraint, find the maximum reactive power from the generator
Q_gen_max_pu = P_gen_pu * tan_phi_max

# The total required reactive power Q_W must be met by the generator and compensator
# Q_W = Q_gen + Q_comp
# To minimize Q_comp, we must maximize Q_gen, so we use Q_gen_max
Q_opt_pu = Q_W_pu - Q_gen_max_pu

# --- Step 7: Final Result ---
# Convert the optimal reactive power from p.u. to MVAR
Q_opt_MVAR = Q_opt_pu * S_base_MVA

print("--- Problem Formulation and Solution ---")
print(f"Base Power (S_base): {S_base_MVA} MVA")
print(f"Target Voltage at Bus-W (V_W): {V_W_pu} p.u.")
print(f"Faulty Line Impedance (Z_fault): {Z_fault_pu:.4f} p.u.")
print("\n--- Power Calculation ---")
print(f"Power required by line to maintain voltage (P_W): {P_W_pu:.4f} p.u.")
print(f"Reactive power required by line (Q_W): {Q_W_pu:.4f} p.u.")
print(f"Total active power from generator including losses (P_gen): {P_gen_pu:.4f} p.u.")
print("\n--- Optimization ---")
print(f"Maximum reactive power from generator at PF={PF_min} (Q_gen_max): {Q_gen_max_pu:.4f} p.u.")
print("\nThe optimal reactive power injection Q_opt is calculated as:")
print("Q_opt = Q_W - Q_gen_max")
print(f"Q_opt (p.u.) = {Q_W_pu:.4f} - {Q_gen_max_pu:.4f}")
print(f"Q_opt (p.u.) = {Q_opt_pu:.4f} p.u.")

print("\n--- Final Answer ---")
print(f"The optimal reactive power injection required from the compensating device is:")
print(f"Q_opt = {Q_opt_MVAR:.2f} MVAR")

# Note on feasibility
if Q_opt_MVAR > Q_max_MVAR:
    print(f"\nNote: The calculated optimal injection ({Q_opt_MVAR:.2f} MVAR) exceeds the device's maximum capacity of {Q_max_MVAR} MVAR.")
    print("This means it's not possible to restore the voltage to the target level while satisfying all constraints with the given hardware.")
else:
    print(f"\nThis value is within the compensator's maximum capacity of {Q_max_MVAR} MVAR.")

print("\nFinal Equation with values:")
print(f"{Q_opt_MVAR:.2f} = ({Q_W_pu:.4f} - {Q_gen_max_pu:.4f}) * {S_base_MVA}")
final_answer = Q_opt_MVAR
# The final answer is the calculated optimal value
# <<<26.43>>>