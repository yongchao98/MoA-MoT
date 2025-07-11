import math

def solve_wind_park_optimization():
    """
    Analyzes the E-STATCOM optimization problem by checking the feasibility of each provided answer choice.
    """
    # 1. System Parameters
    S_base = 10.0  # MVA
    P_wp_MW = 8.0
    P_load_MW = 6.0
    Q_load_MVAR = 2.0
    
    # Line impedance in per unit
    R_pu = 0.05
    X_pu = 0.2
    
    # E-STATCOM limits
    P_ES_max_MW = 4.0
    Q_ES_max_MVAR = 3.0
    S_ES_max_MVA = 5.0
    
    # System Constraints
    V_pcc_min_pu = 0.985
    V_pcc_max_pu = 1.015
    PF_min = 0.98
    
    # Derived constraint for power factor: |Q_g/P_g| <= tan(acos(PF_min))
    tan_phi_max = math.tan(math.acos(PF_min))

    # Options to check
    # Format: [P_ES_MW, Q_ES_MVAR, Stated_Loss_MW]
    options = {
        'A': [1.5, 1.8, 0.5],
        'B': [3.2, 2.1, 0.4],
        'C': [3.0, 2.0, 0.45],
        'D': [2.5, 2.8, 0.35],
        'E': [3.5, 2.5, 0.5]
    }
    
    print("Analyzing each option against the system constraints...")
    print("-" * 60)

    final_choice = None
    
    for key, values in options.items():
        P_ES_MW, Q_ES_MVAR, stated_loss_MW = values
        is_feasible = True
        
        print(f"Option {key}: P_ES = {P_ES_MW} MW, Q_ES = {Q_ES_MVAR} MVAR")

        # Convert to per unit
        P_ES_pu = P_ES_MW / S_base
        Q_ES_pu = Q_ES_MVAR / S_base
        
        # --- Check E-STATCOM Capacity Constraints ---
        s_es_pu = math.sqrt(P_ES_pu**2 + Q_ES_pu**2)
        s_es_max_pu = S_ES_max_MVA / S_base
        
        if abs(P_ES_MW) > P_ES_max_MW or abs(Q_ES_MVAR) > Q_ES_max_MVAR or s_es_pu > s_es_max_pu:
            is_feasible = False
            print(f"  [FAIL] STATCOM capacity exceeded. S_ES={s_es_pu*S_base:.2f} MVA (Limit: {S_ES_max_MVA} MVA)")

        # --- Calculate Grid Power Flow ---
        # P_g = P_wp + P_ES - P_load
        P_g_MW = P_wp_MW + P_ES_MW - P_load_MW
        # Q_g = Q_wp (assumed 0) + Q_ES - Q_load
        Q_g_MVAR = Q_ES_MVAR - Q_load_MVAR
        
        P_g_pu = P_g_MW / S_base
        Q_g_pu = Q_g_MVAR / S_base

        # --- Check Power Factor Constraint ---
        # To avoid division by zero if P_g_pu is 0
        if P_g_pu != 0:
            pf_ratio = abs(Q_g_pu / P_g_pu)
            if pf_ratio > tan_phi_max:
                is_feasible = False
                actual_pf = math.cos(math.atan(pf_ratio))
                print(f"  [FAIL] Power factor constraint not met. Actual PF = {actual_pf:.4f} (Required: > {PF_min})")

        # --- Check Voltage Constraint ---
        # V_pcc ≈ V_grid + R_pu*P_g_pu + X_pu*Q_g_pu  (assuming V_grid = 1.0 pu)
        V_pcc_pu = 1.0 + R_pu * P_g_pu + X_pu * Q_g_pu
        if not (V_pcc_min_pu <= V_pcc_pu <= V_pcc_max_pu):
            is_feasible = False
            print(f"  [FAIL] Voltage constraint not met. V_pcc = {V_pcc_pu:.4f} pu (Range: [{V_pcc_min_pu}, {V_pcc_max_pu}])")
            
        if is_feasible:
            print("  [PASS] All constraints are satisfied.")
            final_choice = key
            # Calculate transmission loss
            # Loss_pu = R_pu * (P_g^2 + Q_g^2) / V_pcc^2. Using simplified Loss_pu = R_pu * (P_g^2 + Q_g^2) is common.
            loss_pu = R_pu * (P_g_pu**2 + Q_g_pu**2)
            calculated_loss_MW = loss_pu * S_base
            
            print(f"\n  Feasible Solution Found: Option {key}")
            print(f"  E-STATCOM Output: P_ES = {P_ES_MW:.2f} MW, Q_ES = {Q_ES_MVAR:.2f} MVAR")
            print(f"  Net Power to Grid: P_g = {P_g_MW:.2f} MW, Q_g = {Q_g_MVAR:.2f} MVAR")
            print(f"  Calculated Transmission Loss = {calculated_loss_MW:.4f} MW")
            print(f"  Stated Total Loss in Option = {stated_loss_MW:.4f} MW")
            print(f"  Note: The calculated transmission loss ({calculated_loss_MW:.4f} MW) does not match the stated loss ({stated_loss_MW} MW), likely due to a typo in the problem's option value or inclusion of unstated E-STATCOM losses.")

        print("-" * 60)

    if final_choice:
        print(f"\nConclusion: Option {final_choice} is the only feasible solution that satisfies all operational constraints.")
    else:
        print("\nConclusion: None of the options are feasible based on the problem statement and standard power system analysis.")
        
solve_wind_park_optimization()
# Final Answer Selection
# Based on the code's output, only Option A satisfies all the constraints.
# Therefore, despite the discrepancy in the loss value, it is the only correct answer among the choices.
final_answer = 'A'
print(f"\nThe final answer is $\\boxed{A}$")