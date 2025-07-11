import cmath

# --- Given Parameters ---
V_actual_pu = 0.85  # Actual voltage at Bus B during fault (pu)
V_target_pu = 1.0   # Target voltage at Bus B (pu)
Rs_pu = 0.02        # System resistance (pu)
Xs_pu = 0.10        # System reactance (pu)
Q_max_mvar = 50.0   # Maximum STATCOM reactive power (MVAR)
S_base_mva = 100.0  # Base power (MVA)
harmonic_loss_increase = 0.04  # 4% increase in losses

# --- Calculations ---

# 1. Convert Q_max to per-unit
Q_max_pu = Q_max_mvar / S_base_mva

# 2. As per the plan, we assume the required reactive power is at the limit of the STATCOM.
# This is a reasoned assumption to solve the under-determined problem.
Q_opt_pu = Q_max_pu

# 3. Calculate the change in voltage and average voltage
delta_V_pu = V_target_pu - V_actual_pu
V_avg_pu = (V_target_pu + V_actual_pu) / 2

# 4. Estimate the total system Thevenin reactance based on the simplified formula
# Q_opt = V_avg * delta_V / X_th  => X_th = V_avg * delta_V / Q_opt
X_th_pu = (V_avg_pu * delta_V_pu) / Q_opt_pu

# 5. Calculate Power Losses
# Assume the Thevenin resistance is approximately Rs_pu, as transformer resistance is small.
R_th_pu = Rs_pu
# Current injected by STATCOM. Assuming purely reactive injection. I = Q/V.
I_grid_pu = Q_opt_pu / V_target_pu
# Base power loss
P_loss_base_pu = R_th_pu * (I_grid_pu ** 2)
# Total power loss including harmonic effects
P_loss_total_pu = P_loss_base_pu * (1 + harmonic_loss_increase)

# 6. Convert results from per-unit to actual values
Q_opt_mvar = Q_opt_pu * S_base_mva
P_loss_total_mw = P_loss_total_pu * S_base_mva

# --- Output the Results ---
print("--- Optimization Results ---")
print(f"To restore the voltage at Bus B to {V_target_pu*100}%, the MMCC STATCOM must inject reactive power.")
print("\nObjective Function: Minimize Q_comp")
print(f"Optimal Reactive Power (Q_opt): {Q_opt_pu:.4f} pu = {Q_opt_mvar:.2f} MVAR")
print(f"This is based on the assumption that the required power equals the STATCOM's maximum capacity.")

print("\nSystem Losses Calculation:")
equation_I = f"I = Q_opt / V_target = {Q_opt_pu:.2f} / {V_target_pu:.2f} = {I_grid_pu:.2f} pu"
print(f"Current from grid: {equation_I}")
equation_Ploss_base = f"P_loss_base = R_th * I^2 = {R_th_pu:.2f} * {I_grid_pu:.2f}^2 = {P_loss_base_pu:.4f} pu"
print(f"Base system losses: {equation_Ploss_base}")
equation_Ploss_total = f"P_loss_total = P_loss_base * (1 + harmonic_factor) = {P_loss_base_pu:.4f} * (1 + {harmonic_loss_increase:.2f}) = {P_loss_total_pu:.4f} pu"
print(f"Total system losses (incl. harmonics): {equation_Ploss_total}")
print(f"Total Real Power Loss: {P_loss_total_pu:.4f} pu = {P_loss_total_mw:.2f} MW")

# Final answer in specified format
# The question asks for Q_opt and the real power losses.
# We will output Q_opt in MVAR.
# <<<50.0>>>
# We will output P_loss in MW.
# <<<0.52>>>