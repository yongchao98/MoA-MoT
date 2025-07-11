import cmath
import math

# System Parameters
S_base = 10.0  # MVA
V_base_pcc = 11.0 # kV
V_base_grid = 220.0 # kV
P_wp_MW = 8.0
Q_wp_MVAR = 0.0 # Assuming wind park operates at unity power factor
P_load_MW = 6.0
Q_load_MVAR = 2.0
Z_line_R_pu = 0.05 # p.u.
Z_line_X_pu = 0.2  # p.u.

# Constraints
P_ES_max_pu = 4.0 / S_base
Q_ES_max_pu = 3.0 / S_base
PF_min = 0.98
V_pcc_min_pu = 0.985
V_pcc_max_pu = 1.015
V_grid_pu = 1.0

# Calculate tan(phi) from PF
tan_phi_max = math.tan(math.acos(PF_min))

# Answer Choices (in MW and MVAR)
choices = {
    'A': {'P_ES': 1.5, 'Q_ES': 1.8, 'Loss': 0.5},
    'B': {'P_ES': 3.2, 'Q_ES': 2.1, 'Loss': 0.4},
    'C': {'P_ES': 3.0, 'Q_ES': 2.0, 'Loss': 0.45},
    'D': {'P_ES': 2.5, 'Q_ES': 2.8, 'Loss': 0.35},
    'E': {'P_ES': 3.5, 'Q_ES': 2.5, 'Loss': 0.5}
}

results = {}

print("Evaluating each option:\n")

for key, values in choices.items():
    P_ES_pu = values['P_ES'] / S_base
    Q_ES_pu = values['Q_ES'] / S_base
    
    # Calculate total injected power at PCC
    P_g_pu = (P_wp_MW / S_base) + P_ES_pu
    Q_g_pu = (Q_wp_MVAR / S_base) + Q_ES_pu
    
    # --- Check Constraints ---
    # 1. Power limits
    p_limit_ok = abs(P_ES_pu) <= P_ES_max_pu
    q_limit_ok = abs(Q_ES_pu) <= Q_ES_max_pu
    
    # 2. Power factor
    if P_g_pu == 0:
        pf_ok = False
    else:
        pf = abs(P_g_pu) / math.sqrt(P_g_pu**2 + Q_g_pu**2)
        pf_ok = pf >= PF_min

    # 3. Voltage
    # Solve v_sq^2 - (V_g^2 + 2*(R*P_g + X*Q_g))*v_sq + |Z|^2*|S|^2 = 0
    # where v_sq = |V_pcc|^2
    R = Z_line_R_pu
    X = Z_line_X_pu
    Z_sq = R**2 + X**2
    S_g_sq = P_g_pu**2 + Q_g_pu**2
    
    A = R * P_g_pu + X * Q_g_pu
    
    # Quadratic equation coeffs: a*x^2 + b*x + c = 0
    a_quad = 1
    b_quad = -(V_grid_pu**2 + 2 * A)
    c_quad = Z_sq * S_g_sq
    
    discriminant = b_quad**2 - 4 * a_quad * c_quad
    v_sq = 0
    if discriminant >= 0:
        # Take the higher voltage solution
        v_sq = (-b_quad + math.sqrt(discriminant)) / (2 * a_quad)
        
    V_pcc_pu = math.sqrt(v_sq) if v_sq > 0 else 0
    v_limit_ok = V_pcc_min_pu <= V_pcc_pu <= V_pcc_max_pu
    
    # --- Calculate Loss ---
    loss_pu = 0
    if v_sq > 0:
        loss_pu = R * S_g_sq / v_sq
    loss_mw = loss_pu * S_base
    
    results[key] = {
        'P_g_pu': P_g_pu,
        'Q_g_pu': Q_g_pu,
        'p_limit_ok': p_limit_ok,
        'q_limit_ok': q_limit_ok,
        'pf_ok': pf_ok,
        'V_pcc_pu': V_pcc_pu,
        'v_limit_ok': v_limit_ok,
        'loss_mw': loss_mw
    }

    # Print evaluation for each choice
    print(f"--- Option {key} ---")
    print(f"P_ES = {values['P_ES']} MW, Q_ES = {values['Q_ES']} MVAR")
    print(f"Calculated Injected Power: P_g = {P_g_pu * S_base:.2f} MW, Q_g = {Q_g_pu * S_base:.2f} MVAR")
    print(f"Constraints Check:")
    print(f"  Power Limits: P {'OK' if p_limit_ok else 'FAIL'}, Q {'OK' if q_limit_ok else 'FAIL'}")
    print(f"  Power Factor: {'OK' if pf_ok else 'FAIL'} (Calculated PF: {pf:.4f})")
    print(f"  Voltage: {'OK' if v_limit_ok else 'FAIL'} (Calculated V_pcc: {V_pcc_pu:.4f} p.u.)")
    print(f"Calculated Loss: {loss_mw:.4f} MW (Stated loss in option: {values['Loss']} MW)\n")


print("--- Conclusion ---")
print("As shown, options D and E fail the power factor constraint.")
print("Options A, B, and C fail the voltage constraint (voltage is too high).")
print("This indicates an inconsistency in the problem's parameters.")
print("Assuming the voltage constraint is relaxed, we should choose the option with the minimum calculated loss among those that satisfy the other constraints (A, B, C).")

min_loss = float('inf')
best_option = None
for key in ['A', 'B', 'C']:
    if results[key]['loss_mw'] < min_loss:
        min_loss = results[key]['loss_mw']
        best_option = key

print(f"\nComparing calculated losses for A, B, C:")
print(f"  Loss(A) = {results['A']['loss_mw']:.3f} MW")
print(f"  Loss(B) = {results['B']['loss_mw']:.3f} MW")
print(f"  Loss(C) = {results['C']['loss_mw']:.3f} MW")
print(f"The minimum loss is calculated for option {best_option}.")
print(f"Therefore, Option {best_option} is the most plausible answer despite the inconsistencies.")

final_choice = choices[best_option]
print(f"\nFinal Answer based on minimizing loss among plausible options:")
print(f"The optimized real and reactive power output of the E-STATCOM are P_ES = {final_choice['P_ES']} MW and Q_ES = {final_choice['Q_ES']} MVAR, with a total system power loss of approximately {results[best_option]['loss_mw']:.2f} MW.")
print(f"This corresponds to option A, which states a loss of {final_choice['Loss']} MW.")

print("\nFinal Equation for Loss Calculation (Option A):")
print(f"P_g = (P_wp + P_ES) = ({P_wp_MW} + {choices['A']['P_ES']}) = {P_wp_MW + choices['A']['P_ES']} MW = {results['A']['P_g_pu']:.2f} p.u.")
print(f"Q_g = (Q_wp + Q_ES) = ({Q_wp_MVAR} + {choices['A']['Q_ES']}) = {Q_wp_MVAR + choices['A']['Q_ES']} MVAR = {results['A']['Q_g_pu']:.2f} p.u.")
print(f"|V_pcc| was calculated to be {results['A']['V_pcc_pu']:.3f} p.u.")
print(f"P_loss = R_pu * (P_g_pu^2 + Q_g_pu^2) / |V_pcc|^2 * S_base")
print(f"P_loss = {Z_line_R_pu} * ({results['A']['P_g_pu']:.2f}^2 + {results['A']['Q_g_pu']:.2f}^2) / {results['A']['V_pcc_pu']:.3f}^2 * {S_base} = {results['A']['loss_mw']:.3f} MW")

<<<A>>>