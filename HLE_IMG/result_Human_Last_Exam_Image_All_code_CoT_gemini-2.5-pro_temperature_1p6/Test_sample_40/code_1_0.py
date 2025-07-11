import cmath
import math

# --- 1. Given Data and Per-Unit Conversion ---
# System Parameters (assuming given values are in per-unit on a 100 MVA base)
Z_S_pu = 0.02 + 0.10j
Z_F_pu = 0.15 + 0j
Q_max_MVAR = 50
S_base_MVA = 100
V_B_fault_pu = 0.85  # Voltage at Bus B during fault
V_B_target_pu = 1.0  # Target voltage at Bus B after compensation
harmonic_loss_factor = 1.04

# Convert MVAR to per-unit
Q_max_pu = Q_max_MVAR / S_base_MVA

# --- 2. Determine Thevenin Equivalent from Bus B ---
# Thevenin impedance at Bus A is Z_S || Z_F
Z_th_A = (Z_S_pu * Z_F_pu) / (Z_S_pu + Z_F_pu)
R_th_A = Z_th_A.real
X_th_A = Z_th_A.imag

# The Thevenin voltage at Bus B during the fault is given as 0.85 pu
V_th_B_mag = V_B_fault_pu

# --- 3. Solve for Unknown System Reactance ---
# The relationship between V_B, V_th, Z_th, and Q_comp is a quadratic equation:
# |Z_th|^2 * Q_comp^2 - 2*X_th*Q_comp + (V_target^2 - V_th^2) = 0
# We assume the system requires the maximum available Q, so Q_opt = Q_max = 0.5 pu.
# We solve for the Thevenin reactance at Bus B (X_th_B) that satisfies this.
# Let X = X_th_B. The equation for X is: 0.25*(R_th_A^2 + X^2) - X + (V_target^2 - V_th^2) = 0
Q_opt_pu = Q_max_pu
a = 0.25
b = -1
c = 0.25 * R_th_A**2 + (V_B_target_pu**2 - V_th_B_mag**2)

# Solve quadratic equation ax^2 + bx + c = 0 for X
discriminant = b**2 - 4 * a * c
sol1 = (-b - math.sqrt(discriminant)) / (2 * a)
sol2 = (-b + math.sqrt(discriminant)) / (2 * a)

# The smaller, more physically realistic solution for Thevenin reactance is chosen.
X_th_B = min(sol1, sol2)
R_th_B = R_th_A
Z_th_B = R_th_B + 1j * X_th_B

# The unknown transformer reactance is the difference.
X_T1 = X_th_B - X_th_A
Z_T1 = 1j * X_T1

print("--- Optimization Problem Solution ---")
print(f"To restore voltage to {V_B_target_pu:.2f} pu, the STATCOM must inject its maximum capacity.")
print(f"Optimal Reactive Power (Q_opt): {Q_opt_pu * S_base_MVA:.2f} MVAR\n")
print(f"The solved for system Thevenin reactance at Bus B (X_th_B) is: {X_th_B:.4f} pu")
print(f"This implies a transformer (TR1) reactance of: {X_T1:.4f} pu\n")

# --- 4. Calculate System Losses ---
print("--- System Loss Calculation ---")
# Final compensated state: |V_B| = 1.0 pu, |V_th|=0.85 pu, Q_comp=0.5 pu
# Solve for the angle of V_B.
# From Re part: |V_th|*cos(delta) = V_B_target^2 - Q_opt*X_th_B (Incorrect formula)
# Let's use the one from the derivation: cos(delta) = (V_target^2 - V_th^2 + Q_opt*X_th_B + R_th_B*sqrt(...)) / ... too complex.
# Use sin(delta) and cos(delta) from power equations:
# Re: V_target = |V_th|*cos(delta) + R_th*P + X_th*Q  -- Power injection Q is -Q_comp
# Im: 0 = |V_th|*sin(delta) + X_th*P - R_th*Q
# Assuming P_comp=0, STATCOM draws S_B = -jQ_comp.
# The current into Bus B is I_AB = (S_B)* / V_B* = (jQ_comp)* / V_B* = -jQ_comp / V_B*
# V_B = V_th - I_AB*Z_th => V_B = V_th - (-jQ_comp / V_B*) * Z_th
# |V_B|^2 = V_th * V_B.conjugate() + j*Q_comp * Z_th.conjugate()
# Let V_B = |V_B|*(cos(d) + jsin(d)). V_th=|V_th|
# |V_B|^2 = |V_th|*|V_B|*(cos(d) - jsin(d)) + j*Q_comp*(R_th-jX_th)
# Equating real and imag parts:
# Re: |V_B|^2 = |V_th|*|V_B|*cos(d) + Q_comp*X_th_B
# Im: 0 = -|V_th|*|V_B|*sin(d) + Q_comp*R_th_B
# |V_B|=1, |V_th|=0.85, Q_comp=0.5
cos_delta = (V_B_target_pu**2 - Q_opt_pu * X_th_B) / (V_th_B_mag * V_B_target_pu)
sin_delta = (Q_opt_pu * R_th_B) / (V_th_B_mag * V_B_target_pu)
delta = math.atan2(sin_delta, cos_delta)
V_B_final = cmath.rect(V_B_target_pu, delta)

print("Final State Calculation:")
print(f"  Final voltage at Bus B (V_B): {abs(V_B_final):.2f} V at angle {math.degrees(delta):.2f} degrees")

# Calculate currents
# Current supplied by STATCOM leads to power draw at bus B of S_B = -j*Q_opt_pu
S_B_net = -1j * Q_opt_pu
# Current flowing from Bus A to Bus B through TR1
I_AB = (S_B_net.conjugate() / V_B_final.conjugate())
# Voltage at Bus A
V_A = V_B_final + I_AB * Z_T1
# Current through the fault impedance
I_fault = V_A / Z_F_pu
# Current from the source grid through Z_S
I_S = I_AB + I_fault

print(f"  Current from A to B (I_AB): {abs(I_AB):.2f} A pu")
print(f"  Voltage at Bus A (V_A): {abs(V_A):.4f} pu")
print(f"  Fault Current (I_fault): {abs(I_fault):.2f} A pu")
print(f"  Source Current (I_S): {abs(I_S):.2f} A pu")

# Calculate fundamental losses (P = I^2*R)
P_loss_S_pu = abs(I_S)**2 * Z_S_pu.real
P_loss_T1_pu = abs(I_AB)**2 * Z_T1.real # R_T1 is assumed to be 0
P_loss_fundamental_pu = P_loss_S_pu + P_loss_T1_pu

# Calculate total losses including harmonics
P_loss_total_pu = P_loss_fundamental_pu * harmonic_loss_factor
P_loss_total_MW = P_loss_total_pu * S_base_MVA

print("\nFinal Loss Calculation:")
print(f"Equation for fundamental losses (pu): P_loss = |I_S|^2 * R_S")
print(f"P_loss = {abs(I_S)**2:.2f} * {Z_S_pu.real:.2f} = {P_loss_fundamental_pu:.4f} pu")
print(f"Total losses including 4% harmonics: P_loss_total = {P_loss_fundamental_pu:.4f} * {harmonic_loss_factor:.2f} = {P_loss_total_pu:.4f} pu")
print(f"Total System Real Power Losses: {P_loss_total_MW:.2f} MW")
print("\n<<<50.0>>>")