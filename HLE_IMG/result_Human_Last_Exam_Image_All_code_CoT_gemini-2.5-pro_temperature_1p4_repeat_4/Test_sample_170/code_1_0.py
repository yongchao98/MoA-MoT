import numpy as np

# Step 1: Define system parameters and state assumptions
S_base = 100.0  # MVA
V_B_nom = 220.0  # kV
V_B_sag_pu = 0.85 # p.u.
V_B_target_pu = 1.0 # p.u.
harmonic_loss_factor = 1.04 # 4% increase

# --- Assumptions ---
# 1. System is represented by a Thevenin equivalent at Bus B.
# 2. Offshore system load is negligible (P_load=0, Q_load=0).
# 3. Thevenin voltage |V_th| is the sagged voltage under no load.
V_th_pu = V_B_sag_pu + 0j
# 4. Thevenin impedance Z_th is assumed to be a plausible value that makes the problem solvable.
R_th_pu = 0.02
X_th_pu = 0.30
# --- End of Assumptions ---

Z_th_pu = R_th_pu + 1j * X_th_pu
Z_th_mag = np.abs(Z_th_pu)
Z_th_angle_rad = np.angle(Z_th_pu)
V_th_mag = np.abs(V_th_pu)

# Step 2 & 3: Solve for Q_opt using power flow equations
# The load at Bus B is S_B = -j*Q_opt. We solve for Q_opt to get |V_B| = 1.0 p.u.

# From the real power equation, solve for the voltage angle difference (delta)
# P_B = 0 => (|V_B|/|V_th|)*cos(beta) = cos(beta - delta)
V_B_mag = V_B_target_pu
cos_beta = np.cos(Z_th_angle_rad)
cos_beta_minus_delta = np.clip((V_B_mag / V_th_mag) * cos_beta, -1.0, 1.0)
beta_minus_delta = np.arccos(cos_beta_minus_delta)

# Choose the physically stable solution (smaller angle deviation)
delta1 = Z_th_angle_rad - beta_minus_delta
delta2 = Z_th_angle_rad + beta_minus_delta
delta = delta1 if np.abs(delta1) < np.abs(delta2) else delta2

# From the reactive power equation, solve for Q_B = -Q_opt
# Q_B = (|V_th|*|V_B|/|Z_th|)*sin(beta - delta) - (|V_B|^2/|Z_th|)*sin(beta)
sin_beta = np.sin(Z_th_angle_rad)
sin_beta_minus_delta = np.sin(Z_th_angle_rad - delta)
Q_B_pu = (V_th_mag * V_B_mag / Z_th_mag) * sin_beta_minus_delta - (V_B_mag**2 / Z_th_mag) * sin_beta
Q_opt_pu = -Q_B_pu
Q_opt_mvar = Q_opt_pu * S_base

# Step 4: Calculate system real power losses
# Loss = R_th * |I_B|^2 where |I_B| = Q_opt_pu / |V_B|
I_B_mag_pu = Q_opt_pu / V_B_mag
P_loss_base_pu = R_th_pu * (I_B_mag_pu**2)

# Step 5: Adjust for harmonic losses
P_loss_total_pu = P_loss_base_pu * harmonic_loss_factor
P_loss_total_mw = P_loss_total_pu * S_base

# Step 6: Print results formatted as equations
print("--- Problem Formulation and Solution ---")
print("System Model (Thevenin Equivalent at Bus B based on assumptions):")
print(f"V_th = {V_th_mag:.2f} p.u.")
print(f"Z_th = {R_th_pu:.2f} + j{X_th_pu:.2f} p.u.")
print("\nSolving Power Flow Equations for Q_opt:")
# Final Equation for Q_opt
print(f"-Q_opt = (|V_th|*|V_B| / |Z_th|) * sin(beta - delta) - (|V_B|^2 / |Z_th|) * sin(beta)")
print(f"-Q_opt = ({V_th_mag:.2f} * {V_B_mag:.2f} / {Z_th_mag:.3f}) * sin({Z_th_angle_rad*180/np.pi:.1f}° - {delta*180/np.pi:.2f}°) - ({V_B_mag:.2f}^2 / {Z_th_mag:.3f}) * sin({Z_th_angle_rad*180/np.pi:.1f}°)")
print(f"-Q_opt = {Q_B_pu:.3f} p.u.")
print(f"Q_opt = {-Q_B_pu:.3f} p.u. = {-Q_B_pu * S_base:.2f} MVAR")

print("\n--- System Real Power Losses Calculation ---")
print(f"P_loss_total = (R_th_pu * Q_opt_pu^2 * S_base) * (1 + 0.04)")
print(f"P_loss_total = ({R_th_pu:.2f} * {-Q_B_pu:.3f}^2 * {S_base:.1f}) * {harmonic_loss_factor:.2f}")
print(f"P_loss_total = {P_loss_total_mw:.2f} MW")
