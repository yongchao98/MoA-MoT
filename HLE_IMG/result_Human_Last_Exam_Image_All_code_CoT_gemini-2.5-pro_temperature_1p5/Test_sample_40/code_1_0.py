import cmath

# --- Givens and Assumptions ---
S_base_MVA = 100.0  # MVA
V_base_kV = 220.0   # kV at Bus B
V_B_fault_pu = 0.85 # p.u. voltage at Bus B during fault
V_B_target_pu = 1.0 # p.u. target voltage
V_th_pu = 1.0       # p.u. Thevenin equivalent voltage source
# Assuming Z_S is given in p.u. as the ohmic value is unrealistic.
# This represents the Thevenin impedance Z_th.
Z_th_pu = complex(0.02, 0.10)
R_th_pu = Z_th_pu.real
X_th_pu = Z_th_pu.imag
Q_max_MVAR = 50.0   # MVAR
harmonic_loss_factor = 0.04

# --- Step 1: Formulate and Solve for K from the fault condition ---
# The voltage at a receiving bus can be calculated with the formula:
# |V_B| = ( |V_th| + sqrt(|V_th|^2 - 4*K) ) / 2
# where K = P_L*R_th + Q_L*X_th
# We solve for K given |V_B| = V_B_fault_pu
print("--- Calculation of required reactive power (Q_opt) ---")
print(f"Step 1: Calculate the effective load-impedance factor K from the voltage sag.")
print(f"The voltage equation is: |V_B| = (|V_th| + sqrt(|V_th|^2 - 4*K)) / 2")
print(f"Given V_B_fault = {V_B_fault_pu} p.u. and V_th = {V_th_pu} p.u., we solve for K:")
# Rearranging the formula for K: K = (|V_th|^2 - (2*|V_B| - |V_th|)^2) / 4
K = (V_th_pu**2 - (2 * V_B_fault_pu - V_th_pu)**2) / 4
print(f"K = ({V_th_pu}^2 - (2*{V_B_fault_pu} - {V_th_pu})^2) / 4 = {K:.4f}")
print("-" * 20)

# --- Step 2: Calculate the required reactive power Q_opt ---
# To restore voltage to V_B_target_pu = 1.0, the new load factor K' must be 0.
# K' = P_L*R_th + (Q_L - Q_opt)*X_th = K - Q_opt*X_th = 0
# This gives the equation for the optimal reactive power Q_opt.
print("Step 2: Calculate the optimal reactive power Q_opt required to restore voltage.")
print(f"To achieve V_B = {V_B_target_pu} p.u., the new factor K' must be 0.")
print(f"The final equation for Q_opt is: Q_opt = K / X_th")
Q_opt_pu = K / X_th_pu
print(f"Substituting the values: Q_opt = {K:.4f} / {X_th_pu:.2f} = {Q_opt_pu:.4f} p.u.")
print("-" * 20)

# --- Step 3: Convert Q_opt to MVAR ---
print("Step 3: Convert Q_opt from p.u. to MVAR.")
Q_opt_MVAR = Q_opt_pu * S_base_MVA
print(f"Q_opt_MVAR = Q_opt_pu * S_base = {Q_opt_pu:.4f} * {S_base_MVA:.1f} = {Q_opt_MVAR:.2f} MVAR")
print(f"\nNote: The required reactive power is {Q_opt_MVAR:.2f} MVAR, which exceeds the STATCOM's maximum capacity of {Q_max_MVAR:.1f} MVAR.")
print("\n" + "="*50 + "\n")

# --- Step 4: Calculate system real power losses ---
print("--- Calculation of system's real power losses ---")
print("Step 4: Determine the system losses in the compensated state.")
# In the compensated state, K' = 0. This implies that the net load at Bus B is zero.
# If the load is zero, the current from the source is zero.
# I_final = (V_th - V_B_target) / Z_th = (1.0 - 1.0) / Z_th = 0
I_final_pu = 0.0
print("For the voltage to be restored to 1.0 p.u., the net load at Bus B must be zero.")
print(f"This means the final current I_final from the source is {I_final_pu} p.u.")

# The fundamental losses are P_loss_fund = |I_final|^2 * R_th
P_loss_fund_pu = (I_final_pu**2) * R_th_pu
P_loss_fund_MW = P_loss_fund_pu * S_base_MVA
print("\nThe fundamental system loss is calculated by: P_loss_fund = |I_final|^2 * R_th")
print(f"P_loss_fund = {I_final_pu}^2 * {R_th_pu:.2f} = {P_loss_fund_pu:.4f} p.u. = {P_loss_fund_MW:.2f} MW")
print("-" * 20)

# --- Step 5: Add harmonic loss effects ---
print("Step 5: Calculate total system losses including harmonic effects.")
print("The final equation for total losses is: P_loss_total = P_loss_fund * (1 + harmonic_loss_factor)")
P_loss_total_MW = P_loss_fund_MW * (1 + harmonic_loss_factor)
print(f"Substituting the values: P_loss_total = {P_loss_fund_MW:.2f} MW * (1 + {harmonic_loss_factor}) = {P_loss_total_MW:.2f} MW")
print("\n" + "="*50)

print("\nFinal Answer:")
print(f"Optimal reactive power required (Q_opt): {Q_opt_MVAR:.2f} MVAR")
print(f"Total system real power losses: {P_loss_total_MW:.2f} MW")