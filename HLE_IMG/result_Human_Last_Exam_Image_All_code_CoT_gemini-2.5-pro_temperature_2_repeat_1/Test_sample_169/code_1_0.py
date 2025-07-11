import cmath
import math

# Step 1: Define system parameters in per-unit (p.u.)
Z_WF_R = 0.01  # p.u.
Z_WF_X = 0.05  # p.u.
R_F = 0.1  # p.u.
V_W_target = 0.575  # p.u. (interpreting "0.575 kV" as a p.u. value)
V_G = 1.0  # p.u. magnitude
S_base = 100  # MVA
Q_max = 10  # MVAR
Q_max_pu = Q_max / S_base # 0.1 p.u.
harmonic_loss_factor = 1.06

print("--- Step 1: System Parameters (p.u.) ---")
print(f"Transmission Impedance Z_WF = {Z_WF_R} + j{Z_WF_X} p.u.")
print(f"Fault Resistance R_F = {R_F} p.u.")
print(f"Target Voltage |V_W| = {V_W_target} p.u.")
print(f"Grid Voltage |V_G| = {V_G} p.u.")
print(f"Harmonic Loss Increase = {(harmonic_loss_factor - 1)*100:.0f}%")
print("-" * 20)

# Step 2: Adjust for harmonic losses
# Resistance is increased to model higher losses.
# Fault resistance is decreased to model higher fault losses (P_loss = V^2/R).
R_line_eff = Z_WF_R * harmonic_loss_factor
R_F_eff = R_F / harmonic_loss_factor
Z_line_eff = complex(R_line_eff, Z_WF_X)

print("\n--- Step 2: Incorporate Harmonic Losses ---")
print(f"Effective Line Resistance R'_line = {R_line_eff:.4f} p.u.")
print(f"Effective Fault Resistance R'_F = {R_F_eff:.4f} p.u.")

# Step 3: Calculate system admittances
Y_line_eff = 1 / Z_line_eff
G_line_eff = Y_line_eff.real
B_line_eff = Y_line_eff.imag

G_F_eff = 1 / R_F_eff

print("\n--- Step 3: Calculate Admittances ---")
print(f"Effective Line Admittance Y'_line = {G_line_eff:.4f} + j({B_line_eff:.4f}) p.u.")
print(f"Effective Fault Conductance G'_F = {G_F_eff:.4f} p.u.")

# Step 4: Solve for the voltage angle delta (δ)
# The real power balance equation is:
# |V_W|^2*(G'_line + G'_F) - |V_W|*|V_G|*(G'_line*cos(δ) + B'_line*sin(δ)) = 0
# Rearranging gives: C = A*cos(δ) + B*sin(δ)
A = V_G * G_line_eff
B = V_G * B_line_eff
C = V_W_target * (G_line_eff + G_F_eff)

# Solve using the R-alpha method: C = R*cos(δ-α)
R = math.sqrt(A**2 + B**2)
alpha = math.atan2(B, A)

cos_delta_minus_alpha = C / R

print("\n--- Step 4: Solve for Angle (δ) ---")
print(f"Solving equation: {C:.4f} = {A:.4f}*cos(δ) + {B:.4f}*sin(δ)")

if abs(cos_delta_minus_alpha) > 1:
    print("Error: No real solution exists for delta. Voltage collapse indicated.")
else:
    delta_minus_alpha1 = math.acos(cos_delta_minus_alpha)
    delta_minus_alpha2 = -math.acos(cos_delta_minus_alpha)

    delta1_rad = delta_minus_alpha1 + alpha
    delta2_rad = delta_minus_alpha2 + alpha
    
    # Choose the angle with the smaller magnitude (more stable solution)
    delta_rad = delta1_rad if abs(delta1_rad) < abs(delta2_rad) else delta2_rad
    delta_deg = math.degrees(delta_rad)
    
    print(f"Found two possible angles for δ: {math.degrees(delta1_rad):.2f}° and {math.degrees(delta2_rad):.2f}°")
    print(f"Choosing the more stable angle: δ = {delta_deg:.2f}°")
    
    # Step 5: Solve for the optimal reactive power Q_opt
    # The imaginary power balance equation is:
    # Q_opt = -|V_W|^2*B'_line - |V_W|*|V_G|*(G'_line*sin(δ) - B'_line*cos(δ))
    
    cos_delta = math.cos(delta_rad)
    sin_delta = math.sin(delta_rad)

    term1 = -V_W_target**2 * B_line_eff
    term2 = -V_W_target * V_G * (G_line_eff * sin_delta - B_line_eff * cos_delta)
    Q_opt_pu = term1 + term2
    Q_opt_MVAR = Q_opt_pu * S_base

    print("\n--- Step 5: Solve for Optimal Reactive Power (Q_opt) ---")
    print("Using the formula: Q_opt = -|V_W|²*B' - |V_W|*|V_G|*(G'*sin(δ) - B'*cos(δ))")
    print("\nFinal Equation with values:")
    
    # Printing the full equation with numbers plugged in
    equation_str = (
        f"Q_opt = -({V_W_target:.3f}**2) * ({B_line_eff:.4f}) - "
        f"{V_W_target:.3f} * {V_G:.1f} * (({G_line_eff:.4f}) * sin({delta_deg:.2f}°) - ({B_line_eff:.4f}) * cos({delta_deg:.2f}°))"
    )
    print(equation_str)
    
    # Recalculating with formatted values for display to match the string
    q_str_part1 = f"Q_opt = -({V_W_target**2:.4f}) * ({B_line_eff:.4f}) - "
    q_str_part2 = f"{V_W_target*V_G:.3f} * (({G_line_eff * sin_delta:.4f}) - ({B_line_eff * cos_delta:.4f}))"
    print(q_str_part1 + q_str_part2)

    q_str_part3 = f"Q_opt = ({term1:.4f}) - {V_W_target*V_G:.3f} * ({G_line_eff * sin_delta - B_line_eff * cos_delta:.4f})"
    print(q_str_part3)
    
    q_str_part4 = f"Q_opt = {term1:.4f} + {term2:.4f}"
    print(q_str_part4)
    
    print("\n--- Final Result ---")
    print(f"The optimal reactive power compensation required is: Q_opt = {Q_opt_pu:.4f} p.u.")
    print(f"This corresponds to {Q_opt_MVAR:.2f} MVAR.")
    
    if Q_opt_pu < 0:
        print("This is reactive power absorption (inductive compensation).")
    else:
        print("This is reactive power injection (capacitive compensation).")

    # Final check against device capacity
    if abs(Q_opt_pu) > Q_max_pu:
        print(f"\nWarning: The required reactive power |{Q_opt_pu:.4f}| p.u. exceeds the device maximum capacity of {Q_max_pu} p.u.")
