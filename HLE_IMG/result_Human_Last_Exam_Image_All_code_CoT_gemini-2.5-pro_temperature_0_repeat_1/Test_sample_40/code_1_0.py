import cmath

# Step 1: Define system parameters and state assumptions
S_base = 100.0  # MVA
V_B_nominal = 220.0  # kV
# Assuming the given impedance values are in per-unit (pu)
Z_S_pu = 0.02 + 0.10j  # pu
Z_F_pu = 0.15 + 0.0j   # pu (assuming resistive fault impedance)
Q_max = 50.0  # MVAR
Q_max_pu = Q_max / S_base # pu
harmonic_loss_increase = 0.04  # 4% increase

# Assumption: The voltage sag to 85% is likely a typo for a solvable problem.
# We assume the voltage sags to 95% of its nominal value.
V_B_fault_pu = 0.95  # pu
V_B_target_pu = 1.0  # pu
V_A_pu = 1.0 # pu, source voltage at Bus A

# Step 2: Analyze the pre-compensation (fault) state to find the equivalent load
# We use the relationship derived from the power flow equation: V_A - V_B = I * Z_S
# (V_A - V_B) * V_B_conj = (P_B1 - jQ_B1) * Z_S
# This gives two equations (real and imaginary parts)
# Re(eq): (V_A - V_B_fault)*V_B_fault = P_B1*R_S + Q_B1*X_S (assuming V_A, V_B are real)
# Im(eq): 0 = P_B1*X_S - Q_B1*R_S  => P_B1 = Q_B1 * (R_S / X_S)

R_S = Z_S_pu.real
X_S = Z_S_pu.imag

# From the imaginary part equation, find the relationship between P_B1 and Q_B1
P_Q_ratio = R_S / X_S

# From the real part equation, solve for Q_B1 and then P_B1
# Note: A more accurate form is (V_A_pu - V_B_fault_pu) * V_B_fault_pu = P_B1*R_S + Q_B1*X_S
# This avoids approximations like V^2.
term = (V_A_pu - V_B_fault_pu) * V_B_fault_pu
Q_B1 = term / (P_Q_ratio * R_S + X_S)
P_B1 = Q_B1 * P_Q_ratio

print("--- Pre-Compensation State Analysis ---")
print(f"Equivalent fault load P_B1 = {P_B1:.4f} pu")
print(f"Equivalent fault load Q_B1 = {Q_B1:.4f} pu\n")

# Step 3: Analyze the post-compensation state to find Q_opt
# To restore voltage to 1.0 pu, the condition P_B2*R_S + Q_B2*X_S = 0 must be met (from the same real part equation where V_A=V_B=1.0)
# P_B2 = P_B1 (STATCOM injects no real power for min Q)
# Q_B2 = Q_B1 - Q_opt
# So, P_B1*R_S + (Q_B1 - Q_opt)*X_S = 0
# P_B1*R_S + Q_B1*X_S = Q_opt*X_S
# Q_opt = (P_B1*R_S + Q_B1*X_S) / X_S

Q_opt_pu = (P_B1 * R_S + Q_B1 * X_S) / X_S
Q_opt_MVAR = Q_opt_pu * S_base

print("--- STATCOM Compensation Calculation ---")
print("The optimization objective is to find the minimum Q_comp to make V_B = 1.0 pu.")
print("This is achieved by setting the STATCOM real power injection to zero (P_comp = 0).")
print("The required reactive power Q_opt is found by solving the power flow equation for the compensated state.")
print(f"Final Equation for Q_opt (pu): Q_opt = (P_B1 * R_S + Q_B1 * X_S) / X_S")
print(f"Q_opt = ({P_B1:.4f} * {R_S} + {Q_B1:.4f} * {X_S}) / {X_S}")
print(f"Optimal Reactive Power Injection (Q_opt): {Q_opt_pu:.4f} pu or {Q_opt_MVAR:.2f} MVAR\n")

# Verify if the solution is within the STATCOM's capacity
if Q_opt_pu > Q_max_pu:
    print(f"Warning: Required Q_opt ({Q_opt_MVAR:.2f} MVAR) exceeds STATCOM capacity ({Q_max:.2f} MVAR).")
else:
    print(f"Required Q_opt ({Q_opt_MVAR:.2f} MVAR) is within STATCOM capacity ({Q_max:.2f} MVAR).\n")


# Step 4: Calculate the system's real power losses after compensation
# Power at Bus B after compensation
P_B2 = P_B1
Q_B2 = Q_B1 - Q_opt_pu

# Current flowing from Bus A to Bus B
# I_S2 = S_B2_conj / V_B2_conj. We need V_B2.
# V_B2 has a small angle theta. sin(theta) = P_B2*X_S - Q_B2*R_S
# Since P_B1*R_S + Q_B1*X_S = Q_opt*X_S, we have P_B2*R_S + Q_B2*X_S = 0.
# This means 1-cos(theta) is approx 0, so theta is very small. We can approximate V_B2 as 1.0+0j.
S_B2_pu = P_B2 + 1j * Q_B2
V_B2_pu = V_B_target_pu + 0j # Approximation
I_S2_pu = S_B2_pu.conjugate() / V_B2_pu.conjugate()

# Base power loss in the transmission line
P_loss_base_pu = (abs(I_S2_pu)**2) * R_S
P_loss_base_MW = P_loss_base_pu * S_base

# Total power loss including harmonic effects
P_loss_total_pu = P_loss_base_pu * (1 + harmonic_loss_increase)
P_loss_total_MW = P_loss_total_pu * S_base

print("--- System Loss Calculation ---")
print("Losses are calculated for the compensated state.")
print(f"Final Equation for Total Losses (MW): P_loss = |I_S2|^2 * R_S * (1 + harmonic_increase) * S_base")
print(f"P_loss = |{I_S2_pu.real:.4f} + j{-I_S2_pu.imag:.4f}|^2 * {R_S} * (1 + {harmonic_loss_increase}) * {S_base}")
print(f"Total System Real Power Losses: {P_loss_total_pu:.6f} pu or {P_loss_total_MW:.4f} MW")

# Final Answer Formatting
print("\n<<<")
print(f"Q_opt = {Q_opt_MVAR:.2f} MVAR")
print(f"P_loss = {P_loss_total_MW:.4f} MW")
print(">>>")