import math

def solve_wind_park_optimization():
    """
    Analyzes the E-STATCOM operation by evaluating given answer choices against system constraints.
    """
    # System Parameters
    S_base = 10.0  # MVA
    P_wp = 8.0  # MW
    Q_wp = 0.0  # MVAR (assumed unity power factor for wind park)
    
    # Line Impedance (assuming given values are per-unit)
    R_pu = 0.05
    X_pu = 0.2
    
    # E-STATCOM Limits
    P_ES_max = 4.0  # MW
    S_ES_max = 5.0  # MVA
    Q_ES_max = math.sqrt(S_ES_max**2 - P_ES_max**2) # MVAR

    # Constraints
    pf_min = 0.98
    # tan(acos(pf)) gives the max Q/P ratio
    q_p_ratio_max = math.tan(math.acos(pf_min))
    v_pcc_deviation = 0.015
    v_grid_pu = 1.0

    # Answer Choices [P_ES (MW), Q_ES (MVAR), Stated Loss (MW)]
    choices = {
        'A': [1.5, 1.8, 0.5],
        'B': [3.2, 2.1, 0.4],
        'C': [3.0, 2.0, 0.45],
        'D': [2.5, 2.8, 0.35],
        'E': [3.5, 2.5, 0.5]
    }
    
    print("Evaluating Answer Choices:\n")
    
    results = {}

    for key, values in choices.items():
        P_ES = values[0]
        Q_ES = values[1]
        stated_loss = values[2]
        
        print(f"--- Checking Option {key}: P_ES = {P_ES} MW, Q_ES = {Q_ES} MVAR ---")

        # Convert to per-unit
        P_ES_pu = P_ES / S_base
        Q_ES_pu = Q_ES / S_base

        # Calculate power at PCC
        P_g = P_wp + P_ES
        Q_g = Q_wp + Q_ES
        P_g_pu = P_g / S_base
        Q_g_pu = Q_g / S_base

        # --- Constraint Checks ---
        
        # 1. E-STATCOM S-rating constraint
        S_ES_sq = P_ES**2 + Q_ES**2
        s_rating_ok = S_ES_sq <= S_ES_max**2
        print(f"1. E-STATCOM Rating Check: P_ES^2 + Q_ES^2 = {S_ES_sq:.2f} <= {S_ES_max**2:.2f} MVA^2? {'OK' if s_rating_ok else 'FAIL'}")

        # 2. Power Factor constraint
        # PF > 0.98  is equivalent to |Q_g / P_g| <= tan(acos(0.98))
        if P_g != 0:
            pf_ok = abs(Q_g / P_g) <= q_p_ratio_max
            print(f"2. Power Factor Check: |Q_g/P_g| = {abs(Q_g/P_g):.3f} <= {q_p_ratio_max:.3f}? {'OK' if pf_ok else 'FAIL'}")
        else:
            pf_ok = False # Undefined PF
            print(f"2. Power Factor Check: P_g is zero. FAIL")
        
        # 3. Voltage constraint
        # V_pcc_pu is approx V_grid_pu + R_pu*P_g_pu + X_pu*Q_g_pu
        v_pcc_pu_approx = v_grid_pu + (R_pu * P_g_pu + X_pu * Q_g_pu)
        v_constraint_ok = (v_grid_pu - v_pcc_deviation) <= v_pcc_pu_approx <= (v_grid_pu + v_pcc_deviation)
        print(f"3. Voltage Check (Note: Likely erroneous in problem statement):")
        print(f"   Calculated V_pcc = {v_pcc_pu_approx:.4f} pu. Required range: [{v_grid_pu - v_pcc_deviation}, {v_grid_pu + v_pcc_deviation}]")
        print(f"   Is voltage constraint met? {'OK' if v_constraint_ok else 'FAIL'}")

        # --- Loss Calculation ---
        # P_loss_pu = R_pu * (P_g_pu^2 + Q_g_pu^2) / V_pcc_pu^2. Assume V_pcc_pu approx 1.0
        P_loss_pu = R_pu * (P_g_pu**2 + Q_g_pu**2)
        calculated_loss_mw = P_loss_pu * S_base
        print(f"Calculated Loss = {calculated_loss_mw:.3f} MW (stated loss was {stated_loss} MW)")
        
        # Store results, ignoring the problematic voltage constraint for feasibility
        results[key] = {
            'feasible': s_rating_ok and pf_ok,
            'loss': calculated_loss_mw,
            'P_ES': P_ES,
            'Q_ES': Q_ES
        }
        print("-" * 35 + "\n")

    # Final Analysis
    print("\n--- Final Analysis ---")
    print("The voltage constraint as formulated makes all options with positive P_ES and Q_ES infeasible.")
    print("Ignoring the voltage constraint to find the most plausible answer among the choices.\n")

    feasible_options = {k: v for k, v in results.items() if v['feasible']}
    
    if not feasible_options:
        print("No feasible options found after ignoring the voltage constraint.")
    else:
        print("Feasible options (ignoring voltage constraint):")
        for key, value in feasible_options.items():
            print(f"  Option {key}: Loss = {value['loss']:.3f} MW")

        # Find the feasible option with the minimum loss
        best_option_key = min(feasible_options, key=lambda k: feasible_options[k]['loss'])
        best_option_data = feasible_options[best_option_key]
        
        print(f"\nThe best choice is Option {best_option_key} as it is feasible (sans voltage constraint) and offers the lowest power loss.")
        
        print(f"\nFinal Answer Details:")
        final_P_ES = best_option_data['P_ES']
        final_Q_ES = best_option_data['Q_ES']
        final_loss = best_option_data['loss']
        
        # Final output of the equation for loss
        final_P_g = P_wp + final_P_ES
        final_Q_g = Q_wp + final_Q_ES
        print(f"The optimized real power output of the E-STATCOM is P_ES = {final_P_ES} MW.")
        print(f"The optimized reactive power output of the E-STATCOM is Q_ES = {final_Q_ES} MVAR.")
        print(f"The total system power loss is calculated as:")
        print(f"P_loss = R_pu * ( (P_wp + P_ES)^2 + (Q_wp + Q_ES)^2 ) * S_base / S_base^2")
        print(f"P_loss = {R_pu} * ( ({P_wp} + {final_P_ES})^2 + ({Q_wp} + {final_Q_ES})^2 ) / {S_base}")
        print(f"P_loss = {R_pu} * ( ({final_P_g})^2 + ({final_Q_g})^2 ) / {S_base}")
        print(f"P_loss = {R_pu} * ( {final_P_g**2:.2f} + {final_Q_g**2:.2f} ) / {S_base}")
        print(f"P_loss = {R_pu} * ( {final_P_g**2 + final_Q_g**2:.2f} ) / {S_base} = {final_loss:.3f} MW")
        print(f"This calculated loss of {final_loss:.3f} MW is very close to the {choices[best_option_key][2]} MW stated in option {best_option_key}.")


solve_wind_park_optimization()
<<<A>>>