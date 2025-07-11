import cmath
import math

# Step 1: Define system parameters and establish the per-unit system
S_base_MVA = 100.0  # MVA
V_base_kV = 0.575  # kV
Z_WF_ohm = 0.01 + 0.05j  # Omega
R_F_ohm = 0.1  # Omega
Q_max_MVAR = 10.0  # MVAR
harmonic_loss_increase = 0.06  # 6% increase
min_power_factor = 0.95
V_W_target_pu = 1.0  # Target voltage is the nominal voltage

S_base_VA = S_base_MVA * 1e6
V_base_V = V_base_kV * 1e3

# Calculate base impedance
Z_base_ohm = V_base_V**2 / S_base_VA
print(f"Base Parameters:")
print(f"S_base = {S_base_MVA} MVA")
print(f"V_base = {V_base_kV} kV")
print(f"Z_base = {Z_base_ohm:.6f} Ohms\n")

# Step 2: Convert system impedances to per-unit
Z_WF_pu = Z_WF_ohm / Z_base_ohm
R_F_pu = R_F_ohm / Z_base_ohm
print(f"Per-Unit Impedances:")
print(f"Z_WF = {Z_WF_pu.real:.4f} + j{Z_WF_pu.imag:.4f} p.u.")
print(f"R_F = {R_F_pu.real:.4f} p.u.\n")

# Step 3 & 4: Model the system and calculate the power demand at the target voltage
# The equivalent impedance of the fault path seen from Bus-W
Z_eq_pu = Z_WF_pu + R_F_pu

# Calculate the complex power consumed by the fault path (S = |V|^2 / Z*)
S_load_pu = V_W_target_pu**2 / Z_eq_pu.conjugate()
P_fund_pu = S_load_pu.real
Q_load_pu = S_load_pu.imag

print("Power Demand Calculation:")
print(f"Equivalent fault path impedance Z_eq = {Z_eq_pu.real:.4f} + j{Z_eq_pu.imag:.4f} p.u.")
print(f"Fundamental complex power demand S_load = {P_fund_pu:.4f} + j{Q_load_pu:.4f} p.u.\n")

# Step 5: Incorporate harmonic losses to find total real power needed
P_total_pu = P_fund_pu * (1 + harmonic_loss_increase)
Q_total_pu = Q_load_pu  # Reactive power demand is unchanged

print("Total Power Supply Calculation (with harmonic losses):")
print(f"Total Real Power P_total = P_fund * (1 + {harmonic_loss_increase}) = {P_fund_pu:.4f} * {1 + harmonic_loss_increase} = {P_total_pu:.4f} p.u.")
print(f"Total Reactive Power Q_total = {Q_total_pu:.4f} p.u.\n")

# Step 6: Apply power factor constraint to find the generator's max reactive power contribution
# To minimize Q_comp, we must maximize Q_gen.
# The max Q_gen is determined by the minimum power factor of the generator.
# PF = P / sqrt(P^2 + Q^2) => Q = P * tan(acos(PF))
tan_phi_max = math.tan(math.acos(min_power_factor))
P_gen_pu = P_total_pu
Q_gen_max_pu = P_gen_pu * tan_phi_max

print("Generator Contribution:")
print(f"Generator must supply all real power: P_gen = {P_gen_pu:.4f} p.u.")
print(f"Max reactive power from generator at PF={min_power_factor}:")
print(f"Q_gen_max = P_gen * tan(acos({min_power_factor})) = {P_gen_pu:.4f} * {tan_phi_max:.4f} = {Q_gen_max_pu:.4f} p.u.\n")

# Step 7: Determine the optimal reactive power injection from the compensator
# Q_opt = Q_total - Q_gen_max
Q_opt_pu = Q_total_pu - Q_gen_max_pu

print("Optimal Compensator Injection Calculation:")
print("The final equation for the optimal reactive power is:")
print(f"Q_opt (p.u.) = Q_total (p.u.) - Q_gen_max (p.u.)")
print(f"Q_opt (p.u.) = {Q_total_pu:.4f} - {Q_gen_max_pu:.4f} = {Q_opt_pu:.4f} p.u.\n")

# Step 8: Convert the result back to MVAR
Q_opt_MVAR = Q_opt_pu * S_base_MVA

print("Final Answer:")
print(f"The optimal reactive power injection Q_opt is {Q_opt_MVAR:.4f} MVAR.")

# Check if the solution is within the device's capacity
if Q_opt_MVAR <= Q_max_MVAR:
    print(f"This is within the compensating device's maximum capacity of {Q_max_MVAR} MVAR.")
else:
    print(f"Warning: This exceeds the compensating device's maximum capacity of {Q_max_MVAR} MVAR.")
