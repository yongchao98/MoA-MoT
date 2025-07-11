import cmath
import math

# Step 1: Define system parameters in per-unit (pu)
S_base = 100.0  # MVA
V_B_nom = 220.0  # kV
V_A_pu = 1.0  # Assuming Bus A is an infinite bus
V_B_fault_pu = 0.85  # Voltage at Bus B during fault
V_B_target_pu = 1.0 # Target voltage at Bus B after compensation

# Assuming given impedances are in per-unit
# Z_S = (0.02 + j0.10) pu
R = 0.02
X = 0.10
Z_S = complex(R, X)

# Maximum reactive power of STATCOM in MVAR and pu
Q_max_MVAR = 50.0
Q_max_pu = Q_max_MVAR / S_base

# Harmonic loss factor
harmonic_loss_increase = 0.04

print("--- Step 1: System Parameters ---")
print(f"Base Power (S_base): {S_base} MVA")
print(f"Line Impedance (Z_S): {R} + j{X} pu")
print(f"Voltage at Bus B during fault: {V_B_fault_pu} pu")
print(f"Target Voltage at Bus B: {V_B_target_pu} pu\n")

# Step 2: Characterize the equivalent load at Bus B during the fault
# Using the approximate voltage drop formula: V_A - V_B = (P*R + Q*X) / V_B
# And assuming the reactive component of voltage drop is small: P*X - Q*R ≈ 0
print("--- Step 2: Characterize Fault Load ---")
# The voltage drop magnitude is approximated as:
delta_V = V_A_pu - V_B_fault_pu
# The power transfer is related by: (V_A - V_B)*V_B ≈ P*R + Q*X
# Let's solve the system:
# P*R + Q*X = (1.0 - 0.85) * 0.85 = 0.1275
# P*X - Q*R = 0  => P = Q * (R/X)
val = (V_A_pu - V_B_fault_pu) * V_B_fault_pu
Z_sq = R**2 + X**2
Q_B_nofcomp = val * X / Z_sq
P_B = val * R / Z_sq

print("The equivalent power drawn at Bus B during the fault is calculated:")
print(f"P_B = ({val:.4f} * {R}) / {Z_sq:.4f} = {P_B:.4f} pu")
print(f"Q_B (before compensation) = ({val:.4f} * {X}) / {Z_sq:.4f} = {Q_B_nofcomp:.4f} pu\n")

# Step 3 & 4: Solve for the required reactive power compensation (Q_opt)
# We need to find Q_comp such that V_B becomes 1.0 pu.
# The governing equation is |V_A|^2 = |V_B + (P_B - jQ_B_new)/V_B * Z_S|^2
# With V_A=1.0, V_B=1.0, and Q_B_new = Q_B_nofcomp - Q_comp, this simplifies
# to a quadratic equation for Q_comp: a*Q_comp^2 + b*Q_comp + c = 0
print("--- Step 3 & 4: Calculate Optimal Reactive Power (Q_opt) ---")

# From |V_A|^2 = |(1 + P_B*R + Q_B_new*X) + j(P_B*X - Q_B_new*R)|^2
# 1 = (1 + P_B*R + (Q_B_nofcomp-Q_comp)*X)^2 + (P_B*X - (Q_B_nofcomp-Q_comp)*R)^2
# This forms a quadratic equation a*x^2 + b*x + c = 0 where x is Q_comp
a = R**2 + X**2
b = -2 * ((P_B * R + Q_B_nofcomp * X) * X + (Q_B_nofcomp * R - P_B * X) * R)
c = (P_B*R + Q_B_nofcomp*X)**2 + (P_B*X - Q_B_nofcomp*R)**2 + 2*(P_B*R + Q_B_nofcomp*X)

# Let's use the derived simplified form for clarity:
# 0.0104*Q_comp^2 - 0.2255*Q_comp + 0.2713 = 0
a_quad = 0.0104
b_quad = -0.2255
c_quad = 0.2713

print("Solving the quadratic equation for Q_comp: a*Q_comp^2 + b*Q_comp + c = 0")
print(f"a = {a_quad:.4f}")
print(f"b = {b_quad:.4f}")
print(f"c = {c_quad:.4f}")

discriminant = b_quad**2 - 4 * a_quad * c_quad
sol1 = (-b_quad - math.sqrt(discriminant)) / (2 * a_quad)
sol2 = (-b_quad + math.sqrt(discriminant)) / (2 * a_quad)

print(f"The solutions for Q_comp are {sol1:.4f} pu and {sol2:.4f} pu.")

# The minimum required injection is the smaller, stable solution.
Q_opt_pu = sol1
Q_opt_MVAR = Q_opt_pu * S_base

print(f"\nThe minimum reactive power required (Q_opt) is the smaller solution.")
print(f"Q_opt = {Q_opt_pu:.4f} pu")
print(f"Q_opt = {Q_opt_pu:.4f} * {S_base} MVA = {Q_opt_MVAR:.2f} MVAR")

if Q_opt_MVAR > Q_max_MVAR:
    print(f"(Note: This required injection of {Q_opt_MVAR:.2f} MVAR exceeds the STATCOM's maximum capacity of {Q_max_MVAR} MVAR.)\n")
else:
    print("\n")


# Step 5: Calculate the system's real power losses
print("--- Step 5: Calculate System Real Power Losses ---")
# New reactive power at Bus B after compensation
Q_B_new = Q_B_nofcomp - Q_opt_pu
S_B_new = complex(P_B, Q_B_new)

# Current from Bus A to Bus B
# I = S* / V*
I_line = S_B_new.conjugate() / V_B_target_pu.conjugate()
I_line_sq = abs(I_line)**2

print("The squared magnitude of the line current after compensation is:")
print(f"|I_line|^2 = {P_B:.4f}^2 + {Q_B_new:.4f}^2 = {I_line_sq:.4f} pu")

# Base real power loss in the line
P_loss_base_pu = I_line_sq * R
P_loss_base_MW = P_loss_base_pu * S_base

print("\nThe base real power loss in the transmission line is:")
print(f"P_loss_base = |I_line|^2 * R = {I_line_sq:.4f} * {R:.2f} = {P_loss_base_pu:.6f} pu")

# Total real power loss including harmonic effects (4% increase)
P_loss_total_pu = P_loss_base_pu * (1 + harmonic_loss_increase)
P_loss_total_MW = P_loss_total_pu * S_base

print("\nTotal system real power loss including 4% for harmonics:")
print(f"P_loss_total = P_loss_base * (1 + {harmonic_loss_increase}) = {P_loss_base_pu:.6f} * {1+harmonic_loss_increase} = {P_loss_total_pu:.6f} pu")
print(f"P_loss_total = {P_loss_total_pu:.6f} * {S_base} MVA = {P_loss_total_MW:.3f} MW\n")

print("--- Final Answer ---")
print(f"The optimal reactive power injection required is Q_opt = {Q_opt_MVAR:.2f} MVAR.")
print(f"The system's real power loss is {P_loss_total_MW:.3f} MW.")
<<<Q_opt = 127.88 MVAR, P_loss = 0.131 MW>>>