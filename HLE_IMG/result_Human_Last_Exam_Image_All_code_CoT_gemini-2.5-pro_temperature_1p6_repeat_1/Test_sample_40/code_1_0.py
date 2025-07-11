import numpy as np
import cmath

# --- 1. Given System Parameters (in per-unit) ---
S_base_MVA = 100.0
V_A = 1.0        # Source voltage (p.u.)
R_S = 0.02       # System resistance (p.u.)
X_S = 0.10       # System reactance (p.u.)
V_B_fault = 0.85 # Voltage at Bus B during fault (p.u.)
V_B_target = 1.0 # Target voltage at Bus B after compensation (p.u.)
harmonic_loss_factor = 1.04 # 4% increase in losses

Z_S = complex(R_S, X_S)
Z_S2 = np.abs(Z_S)**2 # |Z_S|^2

# --- 2. Determine Critical Load (P_L, Q_L) ---
# We assume the faulted voltage corresponds to the voltage stability limit.
# This gives two equations for P_L and Q_L.
# Eq A from stability limit voltage: V_B^2 = (V_A^2 - 2*(R_S*P_L + X_S*Q_L))/2
# 2*V_B_fault^2 = V_A^2 - 2*R_S*P_L - 2*X_S*Q_L
# => 2*R_S*P_L + 2*X_S*Q_L = V_A^2 - 2*V_B_fault^2
A_P = 2 * R_S
A_Q = 2 * X_S
A_C = V_A**2 - 2 * V_B_fault**2

# Eq B from stability limit condition (discriminant of voltage quadratic is zero):
# (2*(R_S*P_L + X_S*Q_L) - V_A^2)^2 - 4*|Z_S|^2*(P_L^2+Q_L^2) = 0
# Substitute from Eq A: (A_C - V_A^2)^2 - 4*Z_S2*(P_L^2+Q_L^2) = 0
B_C = (A_C - V_A**2)**2 / (4 * Z_S2) # This gives P_L^2 + Q_L^2

# Now solve the system of two equations for P_L and Q_L:
# 1) A_P*P_L + A_Q*Q_L = A_C
# 2) P_L^2 + Q_L^2 = B_C
# From (1), P_L = (A_C - A_Q*Q_L) / A_P. Substitute into (2):
# ((A_C - A_Q*Q_L)/A_P)^2 + Q_L^2 = B_C
# (A_C^2 - 2*A_C*A_Q*Q_L + A_Q^2*Q_L^2)/A_P^2 + Q_L^2 = B_C
# (A_Q^2/A_P^2 + 1)*Q_L^2 + (-2*A_C*A_Q/A_P^2)*Q_L + (A_C^2/A_P^2 - B_C) = 0
quad_a = (A_Q/A_P)**2 + 1
quad_b = -2 * A_C * A_Q / (A_P**2)
quad_c = (A_C/A_P)**2 - B_C
q_l_solutions = np.roots([quad_a, quad_b, quad_c])

# Choose the solution representing the offshore wind farm (injecting real power, so P_L < 0)
solutions = []
for Q_L in q_l_solutions:
    P_L = (A_C - A_Q * Q_L) / A_P
    # Wind farm injects power, so P_L at the load bus is negative
    if P_L < 0:
        solutions.append({'P_L': P_L, 'Q_L': Q_L})
# We expect one physically plausible solution
P_L = solutions[0]['P_L']
Q_L = solutions[0]['Q_L']

# --- 3. Calculate Optimal Reactive Power Injection (Q_opt) ---
# With load (P_L, Q_L) known, solve for net reactive power Q_net required to get V_B = 1.0 p.u.
# The power-voltage equation is: 0 = V_B^4 + (2(R*P+X*Q) - V_A^2)V_B^2 + |Z|^2(P^2+Q^2)
# For V_B = 1.0, this simplifies to a quadratic in Q_net = Q_L - Q_opt:
# |Z_S|^2*Q_net^2 + (2*X_S*V_B_target^2 + 2*|Z_S|^2*Q_L_temp)*Q_net + ... This is complex.
# A simpler form is: 0 = Z_S2*P_L^2 + Z_S2*Q_net^2 + 2*R_S*P_L*V_B_target^2 + 2*X_S*Q_net*V_B_target^2 + V_B_target^4 - (V_A*V_B_target)^2
# Since V_A=V_B_target=1, this is: Z_S2*(P_L^2 + Q_net^2) + 2*R_S*P_L + 2*X_S*Q_net = 0
quad_a_qnet = Z_S2
quad_b_qnet = 2 * X_S
quad_c_qnet = Z_S2 * P_L**2 + 2 * R_S * P_L
q_net_solutions = np.roots([quad_a_qnet, quad_b_qnet, quad_c_qnet])

# The two solutions correspond to the stable (high voltage) and unstable (low voltage) operating points.
# We choose the solution with the smaller magnitude of reactive power flow, which corresponds to the stable point.
Q_net = q_net_solutions[np.abs(q_net_solutions).argmin()] if q_net_solutions[0] != q_net_solutions[1] else q_net_solutions[0]

# Calculate Q_opt
Q_opt = Q_L - Q_net
Q_opt_MVAR = Q_opt * S_base_MVA

# --- 4. Calculate System Real Power Losses ---
# Losses are P_loss = |I_line|^2 * R_S
# We need the current, so we find the angle of V_B first. V_B = V_B_target * exp(j*delta)
# from -V_A*V_B*sin(delta) = X_S*P_L - R_S*Q_net
sin_delta = -(X_S * P_L - R_S * Q_net) / (V_A * V_B_target)
cos_delta = np.sqrt(1 - sin_delta**2)

# V_B vector
V_B_vec = V_B_target * complex(cos_delta, sin_delta)
# Line current I = (V_A - V_B) / Z_S
I_line = (V_A - V_B_vec) / Z_S
# Line loss
P_loss = (np.abs(I_line)**2) * R_S

# Include harmonic losses
P_loss_total = P_loss * harmonic_loss_factor
P_loss_total_MW = P_loss_total * S_base_MVA

# --- 5. Print Results ---
print("--- Optimization Results ---")
print(f"Optimal reactive power injection from STATCOM:")
print(f"Q_opt = {Q_opt:.4f} p.u. = {Q_opt_MVAR:.2f} MVAR")
print("\n--- System Losses Analysis ---")
print("Total system real power losses (including 4% harmonic effects):")
print(f"P_loss_total = {P_loss_total:.4f} p.u. = {P_loss_total_MW:.2f} MW")

# For the final output, display the numbers in the final equation.
print("\nFinal Calculation Summary:")
# Q_opt = Q_L - Q_net
print(f"Q_opt ({Q_opt:.4f}) = Q_L ({Q_L:.4f}) - Q_net ({Q_net:.4f})")
# P_loss_total = P_loss * 1.04 = |I|^2 * R_S * 1.04
print(f"P_loss_total ({P_loss_total_MW:.2f} MW) = |I|^2 * R_s * S_base * 1.04 = ({np.abs(I_line):.2f}^2) * {R_S} * {S_base_MVA} * {harmonic_loss_factor}")
