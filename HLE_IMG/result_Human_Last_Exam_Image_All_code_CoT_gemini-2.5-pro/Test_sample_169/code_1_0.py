import numpy as np

# --- 1. Define System Parameters and Assumptions ---

# Base values
V_base_kV = 0.575  # kV
S_base_MVA = 100.0  # MVA
Z_base_Ohm = V_base_kV**2 / S_base_MVA  # Ohm

# Given impedance values in Ohms
R_WF_Ohm = 0.01
X_WF_Ohm = 0.05

# Convert impedances to per-unit (p.u.)
R_WF_pu = R_WF_Ohm / Z_base_Ohm
X_WF_pu = X_WF_Ohm / Z_base_Ohm

# Model harmonic losses by increasing resistance
harmonic_loss_factor = 1.06
R_eff_pu = R_WF_pu * harmonic_loss_factor
Z_eff_pu = R_eff_pu + 1j * X_WF_pu
Z_eff_mag_sq = np.abs(Z_eff_pu)**2

# Assumed and target values
V_W_target_pu = 1.0  # Target voltage at Bus-W (restored to 0.575 kV nominal)
V_th_pu = 0.7  # Assumed Thevenin voltage of the faulted grid
P_W_pu = 0.05  # Assumed active power from wind farm (5 MW)

# Power factor constraint
min_pf = 0.95
# From PF constraint, Q_W <= P_W * tan(acos(PF))
# To minimize Q_comp, we maximize Q_W from the generator.
Q_W_pu = P_W_pu * np.tan(np.arccos(min_pf))

# Compensator capacity
Q_max_MVAR = 10.0
Q_max_pu = Q_max_MVAR / S_base_MVA

# --- 2. Formulate and Solve the Non-linear Equations ---

# The active power flow equation is:
# P_W = (1/|Z|^2) * [R(V_W^2 - V_W*V_th*cos(d)) + X(V_W*V_th*sin(d))]
# This can be rearranged into A*sin(d) + B*cos(d) = C to solve for delta (d)
A = X_WF_pu * V_W_target_pu * V_th_pu
B = -R_eff_pu * V_W_target_pu * V_th_pu
C = P_W_pu * Z_eff_mag_sq - R_eff_pu * V_W_target_pu**2

# Solve for delta using trigonometry
R_trig = np.sqrt(A**2 + B**2)
# Ensure a solution exists
if np.abs(C / R_trig) > 1:
    print("Error: No real solution for delta exists. Power transfer is not possible under these conditions.")
    delta_rad = np.nan
else:
    phi = np.arctan2(B, A)
    delta_rad = np.arcsin(C / R_trig) - phi

# --- 3. Calculate the Required Reactive Power ---

# The reactive power flow equation is:
# Q_total = Q_W + Q_comp = (1/|Z|^2) * [X(V_W^2 - V_W*V_th*cos(d)) - R(V_W*V_th*sin(d))]
if not np.isnan(delta_rad):
    cos_d = np.cos(delta_rad)
    sin_d = np.sin(delta_rad)
    
    term1 = X_WF_pu * (V_W_target_pu**2 - V_W_target_pu * V_th_pu * cos_d)
    term2 = -R_eff_pu * (V_W_target_pu * V_th_pu * sin_d)
    
    Q_total_pu = (term1 + term2) / Z_eff_mag_sq
    
    # Calculate the required compensation
    Q_comp_pu = Q_total_pu - Q_W_pu
    
    # Check against compensator limits
    if Q_comp_pu > Q_max_pu:
        Q_opt_pu = Q_max_pu
    elif Q_comp_pu < 0:
        Q_opt_pu = 0 # Device cannot absorb reactive power
    else:
        Q_opt_pu = Q_comp_pu
        
    Q_opt_MVAR = Q_opt_pu * S_base_MVA
else:
    Q_opt_MVAR = np.nan
    
# --- 4. Print the Final Result and the Equation ---

print("Formulated Nonlinear Optimization Problem:")
print("Minimize: Q_comp")
print("Subject to:")
print(f"  1. V_W = {V_W_target_pu:.3f} p.u. (Target Voltage)")
print(f"  2. Power Factor at Wind Farm >= {min_pf} lagging")
print(f"  3. 0 <= Q_comp <= {Q_max_MVAR} MVAR")
print("  4. Power Flow Equations:")

# Print the final equation with all numbers plugged in
print("\nSolution of the reactive power equation:")
print(f"Q_W + Q_comp = (1/{Z_eff_mag_sq:.5f}) * [{X_WF_pu:.3f}*({V_W_target_pu:.2f}^2 - {V_W_target_pu:.2f}*{V_th_pu:.2f}*cos({delta_rad:.3f})) - {R_eff_pu:.3f}*({V_W_target_pu:.2f}*{V_th_pu:.2f}*sin({delta_rad:.3f}))]")

# Print the intermediate calculated values
print(f"\nCalculated Total Required Reactive Power (Q_total): {Q_total_pu:.4f} p.u.")
print(f"Wind Farm Reactive Power (Q_W): {Q_W_pu:.4f} p.u.")
print(f"Required Compensator Power (Q_comp): {Q_comp_pu:.4f} p.u.")
print("\n--- Final Answer ---")
print(f"The optimal reactive power injection from the fault compensating device is: {Q_opt_MVAR:.4f} MVAR")
print(f"<<<{Q_opt_MVAR:.4f}>>>")
