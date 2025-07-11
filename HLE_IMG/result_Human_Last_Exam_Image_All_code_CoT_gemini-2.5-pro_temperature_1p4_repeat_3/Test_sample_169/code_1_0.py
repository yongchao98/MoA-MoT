import math

# Step 1: Define system parameters and base values
S_base = 100e6  # VA
V_base = 0.575e3 # V
Z_WF = 0.01 + 0.05j # Ohms
R_F = 0.1 # Ohms
V_target_kv = 0.575 # kV
PF_min = 0.95 # lagging
Q_max_mvar = 10 # MVAR
harmonic_loss_increase_factor = 1.06 # 6% increase

# Calculate base impedance
Z_base = V_base**2 / S_base

# Convert parameters to per-unit (p.u.)
Z_WF_pu = Z_WF / Z_base
R_F_pu = R_F / Z_base
V_target_pu = (V_target_kv * 1000) / V_base
Q_max_pu = (Q_max_mvar * 1e6) / S_base

# Step 2: Model the equivalent impedance of the faulted line
# Z_eq is the line impedance to the fault + fault resistance
Z_eq_pu = Z_WF_pu + R_F_pu

# Step 3: Calculate the complex power drawn by the faulted line
# S = V^2 / Z*  (where Z* is the complex conjugate of Z)
# We need to supply this power to maintain V_target_pu at Bus-W.
S_line_pu = V_target_pu**2 / Z_eq_pu.conjugate()
P_line_pu = S_line_pu.real
Q_line_pu = S_line_pu.imag

# Step 4: Account for harmonic losses
# Total active power required from the generator is increased by 6%
P_gen_pu = P_line_pu * harmonic_loss_increase_factor
# The reactive power requirement of the line is assumed to remain Q_line_pu
Q_needed_pu = Q_line_pu

# Step 5: Apply the power factor constraint to the generator
# PF = cos(phi) >= 0.95 -> |tan(phi)| <= tan(acos(0.95))
# Q_gen / P_gen <= tan(acos(0.95))
# To minimize Q_comp, we must maximize Q_gen
tan_phi_max = math.tan(math.acos(PF_min))
Q_gen_max_pu = P_gen_pu * tan_phi_max

# Step 6: Determine the optimal reactive power injection from the compensator
# Q_needed = Q_gen + Q_comp
# Q_comp_min = Q_needed - Q_gen_max
Q_opt_pu = Q_needed_pu - Q_gen_max_pu

# Step 7: Convert the final answer to MVAR and verify constraints
Q_opt_mvar = Q_opt_pu * S_base / 1e6

# Print the results and the final equation
print("--- System Parameters (in p.u.) ---")
print(f"Target Voltage (V_target): {V_target_pu:.4f} p.u.")
print(f"Equivalent Faulted Line Impedance (Z_eq): ({Z_eq_pu.real:.4f} + {Z_eq_pu.imag:.4f}j) p.u.")
print(f"Compensator Max Capacity (Q_max): {Q_max_pu:.4f} p.u.\n")

print("--- Calculation of Optimal Reactive Power ---")
print("The optimal reactive power injection (Q_opt) is determined by the equation:")
print("Q_opt = Q_needed - Q_gen_max\n")

print("1. Calculate total reactive power needed by the line (Q_needed):")
print(f"   Q_needed = Im(V_target^2 / Z_eq*) = {Q_needed_pu:.5f} p.u.\n")

print("2. Calculate maximum reactive power from the generator (Q_gen_max) under PF constraints:")
print(f"   Active power from generator (P_gen, incl. harmonic losses) = {P_gen_pu:.5f} p.u.")
print(f"   Q_gen_max = P_gen * tan(acos({PF_min})) = {Q_gen_max_pu:.5f} p.u.\n")

print("3. Solve for the optimal reactive power injection (Q_opt):")
print(f"   Q_opt (p.u.) = {Q_needed_pu:.5f} - {Q_gen_max_pu:.5f} = {Q_opt_pu:.5f} p.u.\n")

# Check if the result is within the compensator's capacity
if Q_opt_mvar <= Q_max_mvar:
    print(f"The calculated optimal reactive power injection is {Q_opt_mvar:.3f} MVAR.")
    print(f"This is within the device's maximum capacity of {Q_max_mvar} MVAR.")
else:
    print(f"The required reactive power {Q_opt_mvar:.3f} MVAR exceeds the device's capacity.")

<<<0.264>>>