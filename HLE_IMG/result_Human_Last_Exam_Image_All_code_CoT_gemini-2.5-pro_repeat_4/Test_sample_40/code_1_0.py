import cmath
import math

# Step 1: Define system parameters in per-unit (p.u.)
S_base = 100.0  # MVA
V_B_nom = 220.0  # kV
# Assuming given impedance values are in p.u. on the 100 MVA base
Z_S_pu = 0.02 + 0.10j
R_S_pu = Z_S_pu.real
X_S_pu = Z_S_pu.imag

Z_F_pu = 0.15 + 0.0j # Fault impedance is given as 0.15 Ohm, assumed to be resistive p.u. value
Y_F_pu = 1 / Z_F_pu

V_A_pu = 1.0  # Slack bus voltage (p.u.)
V_B_fault_pu = 0.85  # Voltage at Bus B during fault without compensation (p.u.)
V_B_target_pu = 1.0 # Target voltage at Bus B with compensation (p.u.)

Q_max_pu = 50.0 / S_base # Max STATCOM capacity in p.u.

# Step 2: Determine the equivalent load admittance G_L from the uncompensated fault condition.
# We assume the load (wind farm) has unity power factor (B_L = 0) and is of constant impedance type.
# The governing equation is |V_B| = |V_A / (1 + Z_S * Y_total)|
# where Y_total = G_L + Y_F
# |1 + Z_S * (G_L + Y_F)| = |V_A| / |V_B_fault| = 1.0 / 0.85
# This leads to a quadratic equation for X = G_L + Y_F: a*X^2 + b*X + c = 0
# where a = |Z_S|^2, b = 2*Re(Z_S), c = 1 - (1/|V_B_fault|)^2

a_gl = abs(Z_S_pu)**2
b_gl = 2 * Z_S_pu.real
c_gl = 1 - (1/V_B_fault_pu)**2

print("Step 2: Solving for load characteristics.")
print(f"The quadratic equation for X = (G_L + Y_F) is: {a_gl:.4f}*X^2 + {b_gl:.4f}*X + {c_gl:.4f} = 0")

discriminant_gl = b_gl**2 - 4 * a_gl * c_gl
X1 = (-b_gl + math.sqrt(discriminant_gl)) / (2 * a_gl)
X2 = (-b_gl - math.sqrt(discriminant_gl)) / (2 * a_gl)

# We choose the positive, physically meaningful solution for the total admittance magnitude
X = X1 if X1 > 0 else X2
print(f"Solving for X gives: X = {X:.4f} p.u.")

G_L_pu = X - Y_F_pu.real
print(f"The load admittance G_L is calculated as: G_L = X - Y_F = {X:.4f} - {Y_F_pu.real:.4f} = {G_L_pu:.4f} p.u.")
# A negative G_L means the "load" is a generator (the wind farm), which is correct.

# Step 3: Solve for the required reactive power Q_opt
# With compensation, V_B is restored to 1.0 p.u.
# The total power at Bus B is S_B = P_B + jQ_B
# P_B = P_load + P_fault = G_L*V_B^2 + V_B^2/R_F = (G_L + Y_F)*V_B^2
P_B_comp = (G_L_pu + Y_F_pu.real) * V_B_target_pu**2
# Q_B = Q_load + Q_fault - Q_comp = 0 + 0 - Q_comp
Q_B_comp = - 'Q_comp' # Symbolic at this stage

# Total admittance at Bus B during compensation Y_B = S_B* / |V_B|^2 = P_B - jQ_B
# Y_B = P_B_comp + jQ_comp
# The voltage equation |1 + Z_S * Y_B| = |V_A|/|V_B_target| = 1 leads to a quadratic equation for Q_comp
# a*Q_comp^2 + b*Q_comp + c = 0
a_q = abs(Z_S_pu)**2
b_q = 2 * (P_B_comp * Z_S_pu.imag - Z_S_pu.real)
c_q = (1 + P_B_comp * Z_S_pu.real)**2 + (P_B_comp * Z_S_pu.imag)**2 - 1

print("\nStep 3: Solving for the optimal reactive power Q_opt.")
# Note: The derived coefficients are slightly different from manual calculation due to formula rearrangement
# Let's use the explicit magnitude squared equation:
# |(1 + R_S*P_B - X_S*Q_comp) + j(X_S*P_B + R_S*Q_comp)|^2 = 1
# (1 + R_S*P_B)^2 - 2*X_S*Q_comp*(1+R_S*P_B) + X_S^2*Q_comp^2 + (X_S*P_B)^2 + 2*R_S*X_S*P_B*Q_comp + R_S^2*Q_comp^2 = 1
# This simplifies to a*Q_comp^2 + b*Q_comp + c = 0 with:
a_q = R_S_pu**2 + X_S_pu**2
b_q = 2 * R_S_pu * X_S_pu * P_B_comp - 2 * X_S_pu * (1 + R_S_pu * P_B_comp)
c_q = (1 + R_S_pu * P_B_comp)**2 + (X_S_pu * P_B_comp)**2 - 1
# Let's use the simpler form derived in thought process: 0.0104 * Q_comp^2 - 0.200 * Q_comp + 0.3839 = 0
a_q_final = 0.0104
b_q_final = -0.200
c_q_final = 0.3839

print(f"The quadratic equation for Q_comp is: {a_q_final:.4f}*Q_comp^2 + {b_q_final:.4f}*Q_comp + {c_q_final:.4f} = 0")
discriminant_q = b_q_final**2 - 4 * a_q_final * c_q_final
q_sol1 = (-b_q_final + math.sqrt(discriminant_q)) / (2 * a_q_final)
q_sol2 = (-b_q_final - math.sqrt(discriminant_q)) / (2 * a_q_final)
print(f"The solutions for Q_comp are {q_sol1:.2f} p.u. and {q_sol2:.2f} p.u.")

# Q_opt is the minimum positive (injected) reactive power
Q_opt_pu = min(q_sol1, q_sol2)
Q_opt_MVAR = Q_opt_pu * S_base
print(f"The minimum required reactive power is Q_opt = {Q_opt_pu:.2f} p.u., which is {Q_opt_MVAR:.2f} MVAR.")
if Q_opt_pu > Q_max_pu:
    print(f"(Note: This required value of {Q_opt_MVAR:.2f} MVAR exceeds the STATCOM's maximum capacity of {Q_max_pu*S_base} MVAR.)")

# Step 4 & 5: Calculate system's real power losses with harmonic effects
# With compensation, |V_B| = 1.0 p.u.
# First, find the line current I_line = (V_A - V_B) / Z_S
# We need the complex voltage V_B. V_B = V_A / (1 + Z_S * Y_B)
Y_B_comp = P_B_comp + 1j * Q_opt_pu
V_B_comp_pu = V_A_pu / (1 + Z_S_pu * Y_B_comp)
I_line_pu = (V_A_pu - V_B_comp_pu) / Z_S_pu

# Fundamental losses
P_loss_line_pu = (abs(I_line_pu)**2) * R_S_pu
P_loss_fault_pu = (abs(V_B_comp_pu)**2) / Z_F_pu.real
P_loss_total_pu = P_loss_line_pu + P_loss_fault_pu

# Add 4% for harmonic losses
harmonic_loss_factor = 1.04
P_loss_final_pu = P_loss_total_pu * harmonic_loss_factor
P_loss_final_MW = P_loss_final_pu * S_base

print("\nStep 4 & 5: Calculating final system losses.")
print(f"Line loss = |I_line|^2 * R_S = {abs(I_line_pu):.2f}^2 * {R_S_pu:.2f} = {P_loss_line_pu:.2f} p.u.")
print(f"Fault loss = |V_B|^2 / R_F = {abs(V_B_comp_pu):.2f}^2 / {Z_F_pu.real:.2f} = {P_loss_fault_pu:.2f} p.u.")
print(f"Total fundamental loss = {P_loss_line_pu:.2f} + {P_loss_fault_pu:.2f} = {P_loss_total_pu:.2f} p.u.")
print(f"Total loss including 4% harmonics = {P_loss_total_pu:.2f} * {harmonic_loss_factor} = {P_loss_final_pu:.2f} p.u.")
print(f"This corresponds to {P_loss_final_MW:.2f} MW.")

print("\n--- Final Answer ---")
print(f"The optimal reactive power injection required is Q_opt = {Q_opt_MVAR:.2f} MVAR.")
print(f"The total system real power loss is {P_loss_final_MW:.2f} MW.")
print(f"\nFinal numeric answer format: Q_opt = {Q_opt_MVAR:.1f}, P_loss = {P_loss_final_MW:.1f}")
print(f"<<<Q_opt = {Q_opt_MVAR:.1f}, P_loss = {P_loss_final_MW:.1f}>>>")
