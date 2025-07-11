import math

def evaluate_options():
    # System Parameters
    S_base = 10.0  # MVA
    V_pcc_nominal = 1.0 # pu
    
    # Wind Park
    P_wp = 8.0  # MW
    P_wp_pu = P_wp / S_base
    
    # Line Impedance (assuming given values are in pu)
    R_line_pu = 0.05
    X_line_pu = 0.2
    
    # Constraints
    P_es_max_pu = 4.0 / S_base
    Q_es_max_pu = 3.0 / S_base
    pf_min = 0.98
    v_pcc_min = 0.985 # 1.0 - 0.015
    v_pcc_max = 1.015 # 1.0 + 0.015
    
    # Power factor constraint derived value
    pf_angle_max_rad = math.acos(pf_min)
    q_p_ratio_max = math.tan(pf_angle_max_rad)

    # Answer Choices [P_es_mw, Q_es_mvar, Loss_mw]
    options = {
        'A': [1.5, 1.8, 0.5],
        'B': [3.2, 2.1, 0.4],
        'C': [3.0, 2.0, 0.45],
        'D': [2.5, 2.8, 0.35],
        'E': [3.5, 2.5, 0.5]
    }
    
    results = {}
    
    print("--- Evaluating Options ---")
    
    for key, values in options.items():
        p_es_mw, q_es_mvar, stated_loss_mw = values
        
        # Convert to per-unit
        p_es_pu = p_es_mw / S_base
        q_es_pu = q_es_mvar / S_base
        
        # Calculate power at PCC
        p_g_pu = P_wp_pu + p_es_pu
        q_g_pu = q_es_pu # Assuming Q_wp = 0
        
        # --- Check Constraints ---
        # 1. E-STATCOM power limits
        p_es_feasible = abs(p_es_pu) <= P_es_max_pu
        q_es_feasible = abs(q_es_pu) <= Q_es_max_pu
        
        # 2. Power factor constraint
        # Avoid division by zero if p_g_pu is zero
        pf_feasible = (p_g_pu != 0) and (abs(q_g_pu / p_g_pu) <= q_p_ratio_max)

        # 3. Voltage constraint
        v_pcc_pu_approx = V_pcc_nominal + p_g_pu * R_line_pu + q_g_pu * X_line_pu
        v_feasible = v_pcc_min <= v_pcc_pu_approx <= v_pcc_max
        
        # Calculate objective function (proportional to loss)
        objective_value = p_g_pu**2 + q_g_pu**2
        
        # Calculate line loss based on this option
        # Assuming V_pcc is close to nominal for calculation, as the calculated voltage is out of spec
        calculated_loss_pu = R_line_pu * objective_value / (V_pcc_nominal**2)
        calculated_loss_mw = calculated_loss_pu * S_base
        
        results[key] = {
            "P_ES (pu)": p_es_pu,
            "Q_ES (pu)": q_es_pu,
            "P_g (pu)": p_g_pu,
            "Q_g (pu)": q_g_pu,
            "P_ES limit feasible": p_es_feasible,
            "Q_ES limit feasible": q_es_feasible,
            "PF feasible": pf_feasible,
            "Voltage feasible": v_feasible,
            "V_pcc_approx (pu)": v_pcc_pu_approx,
            "Overall Feasible (excl. Voltage)": p_es_feasible and q_es_feasible and pf_feasible,
            "Objective (P_g^2+Q_g^2)": objective_value,
            "Calculated Loss (MW)": calculated_loss_mw,
            "Stated Loss (MW)": stated_loss_mw
        }

    # Print analysis for each option
    for key, res in results.items():
        print(f"\nOption {key}: P_ES = {options[key][0]} MW, Q_ES = {options[key][1]} MVAR")
        print(f"  P_g = {res['P_g (pu)']*S_base:.2f} MW, Q_g = {res['Q_g (pu)']*S_base:.2f} MVAR")
        print(f"  Constraints Check:")
        print(f"    - P_ES limit: {'OK' if res['P_ES limit feasible'] else 'FAIL'}")
        print(f"    - Q_ES limit: {'OK' if res['Q_ES limit feasible'] else 'FAIL'}")
        print(f"    - Power Factor: {'OK' if res['PF feasible'] else 'FAIL'} (Ratio |Q/P|={abs(res['Q_g (pu)']/res['P_g (pu)']):.3f} vs max {q_p_ratio_max:.3f})")
        print(f"    - Voltage: {'OK' if res['Voltage feasible'] else 'FAIL'} (V_pcc~{res['V_pcc_approx (pu)']:.4f} pu, limit [{v_pcc_min}, {v_pcc_max}])")
        print(f"  Loss Analysis:")
        print(f"    - Objective function value (to be minimized): {res['Objective (P_g^2+Q_g^2)']:.4f}")
        print(f"    - Calculated line loss: {res['Calculated Loss (MW)']:.3f} MW")
        print(f"    - Stated total loss: {res['Stated Loss (MW)']:.3f} MW")

    print("\n--- Conclusion ---")
    print("The voltage constraint is violated by all options, indicating an inconsistency in the problem statement's data.")
    print("Ignoring the voltage constraint, we check the feasibility of the remaining constraints.")
    
    feasible_options = {k: v for k, v in results.items() if v['Overall Feasible (excl. Voltage)']}
    
    if not feasible_options:
        print("No options are feasible based on E-STATCOM and Power Factor limits.")
    else:
        print(f"Feasible options (ignoring voltage): {list(feasible_options.keys())}")
        best_option_key = min(feasible_options, key=lambda k: feasible_options[k]['Objective (P_g^2+Q_g^2)'])
        print(f"The objective is to minimize line losses, which is proportional to P_g^2 + Q_g^2.")
        print(f"Comparing the objective function values for the feasible options:")
        for key, res in feasible_options.items():
            print(f"  - Option {key}: Objective = {res['Objective (P_g^2+Q_g^2)']:.4f}")
        
        print(f"\nOption {best_option_key} has the minimum objective function value, making it the most optimal choice among the given options.")
        final_p = options[best_option_key][0]
        final_q = options[best_option_key][1]
        final_loss = options[best_option_key][2]
        print(f"The final answer is that the optimized real and reactive power output of the E-STATCOM are P_ES = {final_p} MW and Q_ES = {final_q} MVAR, with a total system power loss of {final_loss} MW.")

evaluate_options()
<<<A>>>