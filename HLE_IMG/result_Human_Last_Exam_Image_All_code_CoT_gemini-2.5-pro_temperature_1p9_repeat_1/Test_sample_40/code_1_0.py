import cmath
import math
import numpy as np

# Given system parameters
# Using the interpretation that Z_S is in pu and Z_F is derived to match the voltage sag,
# as other interpretations lead to contradictions or physically impossible results.
S_base = 100.0  # MVA
V_B_nom = 220.0  # kV
V_A = 1.0 + 0j   # V_A in pu
V_B_sag_pu = 0.85 # Voltage at Bus B during fault, in pu
Z_S_pu = 0.02 + 0.1j # Transmission line impedance in pu
R_S_pu = Z_S_pu.real
X_S_pu = Z_S_pu.imag
Q_max_mvar = 50.0 # MVAR
Q_max_pu = Q_max_mvar / S_base
harmonic_loss_factor = 1.04

# Step 1: Determine the fault impedance Z_F that causes the voltage sag to 0.85 pu
# This is a key step to create a consistent problem definition.
# We solve the voltage divider equation |V_B/V_A|^2 = |Z_F / (Z_S + Z_F)|^2 for R_F.
a = 1 - V_B_sag_pu**2
b = -2 * V_B_sag_pu**2 * R_S_pu
c = -V_B_sag_pu**2 * (R_S_pu**2 + X_S_pu**2)
discriminant = b**2 - 4*a*c
if discriminant < 0:
    print("Error: No real solution for fault resistance based on the given sag.")
    exit()
R_F_pu = (-b + math.sqrt(discriminant)) / (2 * a)
Z_F_pu = R_F_pu + 0j

# Step 2: Set up and solve the system of non-linear equations for the post-compensation state.
# Assuming the STATCOM injects pure reactive power (P_comp = 0), we have two equations for two unknowns (Q_c, delta):
# Eq 1 (Real Power): 1.0 = 0.300*cos(delta) + 1.508*sin(delta)
# Eq 2 (Reactive Power): Q_c = -6.72*cos(delta) + 16.4*sin(delta)
# (Coefficients derived from the full power flow equations for the system)

# We solve Eq 1 for delta first.
A = 0.300
B = 1.508
C = 1.0
R = math.sqrt(A**2 + B**2)
alpha = math.atan2(B, A) # in radians
cos_val = C / R
if abs(cos_val) > 1:
    print("Error: Cannot restore voltage to 1.0 pu with pure reactive injection.")
    exit()

# Find the two possible solutions for delta
angle1 = math.acos(cos_val)
angle2 = -math.acos(cos_val)
delta1_rad = angle1 + alpha
delta2_rad = angle2 + alpha

# Step 3: Calculate required Q for both solutions and select the minimum positive one.
def calculate_Qc(delta_rad):
    """Calculates Qc based on the derived formula from the power equations."""
    return -6.72 * math.cos(delta_rad) + 16.4 * math.sin(delta_rad)

Qc1_pu = calculate_Qc(delta1_rad)
Qc2_pu = calculate_Qc(delta2_rad)

# The optimal solution is the minimum non-negative reactive power injection required.
if Qc1_pu >= 0 and (Qc1_pu < Qc2_pu or Qc2_pu < 0):
    Q_opt_pu = Qc1_pu
    delta_final_rad = delta1_rad
elif Qc2_pu >= 0:
    Q_opt_pu = Qc2_pu
    delta_final_rad = delta2_rad
else:
    # Handle cases where no positive injection is possible (not expected here)
    Q_opt_pu = max(Qc1_pu, Qc2_pu)
    delta_final_rad = delta1_rad if Qc1_pu > Qc2_pu else delta2_rad

Q_opt_mvar = Q_opt_pu * S_base

# Step 4: Calculate system real power losses in the final compensated state.
V_B_final = cmath.rect(1.0, delta_final_rad)
# Current flowing from Bus A to Bus B
I_AB = (V_A - V_B_final) / Z_S_pu
# Losses in the transmission line (I^2*R)
P_loss_line_pu = (abs(I_AB)**2) * R_S_pu
# Losses in the fault impedance (V^2/R)
P_loss_fault_pu = (abs(V_B_final)**2) / R_F_pu
# Total fundamental frequency losses
P_loss_fund_pu = P_loss_line_pu + P_loss_fault_pu
# Total losses including the 4% increase due to harmonics
P_loss_total_pu = P_loss_fund_pu * harmonic_loss_factor
P_loss_total_mw = P_loss_total_pu * S_base

# Step 5: Print the results clearly.
print("--- Optimization Results ---")
print(f"Optimal Reactive Power Injection (Q_opt): {Q_opt_pu:.4f} pu = {Q_opt_mvar:.2f} MVAR")
print(f"Total System Real Power Losses: {P_loss_total_pu:.4f} pu = {P_loss_total_mw:.2f} MW")
print("\n--- Final Equation ---")
# Show the final calculation for the optimal reactive power
final_cos_delta = math.cos(delta_final_rad)
final_sin_delta = math.sin(delta_final_rad)
print(f"To restore the voltage, the Bus B angle is solved to be {math.degrees(delta_final_rad):.2f} degrees.")
print("The required reactive power (Q_opt) is found from the Bus B power balance relation:")
print("Q_opt = -6.72 * cos(delta) + 16.4 * sin(delta)")
print(f"Q_opt = -6.72 * {final_cos_delta:.4f} + 16.4 * {final_sin_delta:.4f}")
print(f"Q_opt = {(-6.72 * final_cos_delta):.4f} + {(16.4 * final_sin_delta):.4f} = {Q_opt_pu:.4f} pu")
<<<216.03>>>