import math

def solve_wind_park_optimization():
    """
    Solves the wind park optimization problem by checking each multiple-choice option
    against the given constraints.
    """
    # System Parameters
    P_wp = 8.0  # MW
    Q_wp = 0.0  # MVAR (assuming unity power factor for the wind park)

    # E-STATCOM Constraints
    P_es_max = 4.0  # MW
    Q_es_max = 3.0  # MVAR
    S_es_max = 5.0  # MVA

    # PCC Constraints
    pf_min = 0.98
    # The voltage constraint is assumed to be met due to negligible line impedance,
    # as calculated in the analysis.

    # Answer Choices
    # Format: {'P_ES': MW, 'Q_ES': MVAR, 'Loss': MW}
    options = {
        'A': {'P_ES': 1.5, 'Q_ES': 1.8, 'Loss': 0.5},
        'B': {'P_ES': 3.2, 'Q_ES': 2.1, 'Loss': 0.4},
        'C': {'P_ES': 3.0, 'Q_ES': 2.0, 'Loss': 0.45},
        'D': {'P_ES': 2.5, 'Q_ES': 2.8, 'Loss': 0.35},
        'E': {'P_ES': 3.5, 'Q_ES': 2.5, 'Loss': 0.5},
    }

    valid_options = []
    
    print("Evaluating each option against the constraints...\n")

    for key, values in options.items():
        P_es = values['P_ES']
        Q_es = values['Q_ES']
        loss = values['Loss']
        
        print(f"--- Checking Option {key} ---")
        print(f"P_ES = {P_es} MW, Q_ES = {Q_es} MVAR, Total Loss = {loss} MW")

        # 1. Check E-STATCOM Power Limits
        S_es = math.sqrt(P_es**2 + Q_es**2)
        
        statcom_p_ok = P_es <= P_es_max
        statcom_q_ok = abs(Q_es) <= Q_es_max
        statcom_s_ok = S_es <= S_es_max
        
        statcom_ok = statcom_p_ok and statcom_q_ok and statcom_s_ok

        print(f"E-STATCOM S limit check: {S_es:.2f} MVA <= {S_es_max} MVA -> {'OK' if statcom_s_ok else 'FAIL'}")
        
        # 2. Check PCC Power Factor Constraint
        P_g = P_wp + P_es
        Q_g = Q_wp + Q_es
        S_g = math.sqrt(P_g**2 + Q_g**2)
        
        pf_pcc = P_g / S_g if S_g > 0 else 0
        pf_ok = pf_pcc >= pf_min
        
        print(f"PCC Power Factor check: {pf_pcc:.4f} >= {pf_min} -> {'OK' if pf_ok else 'FAIL'}")

        if statcom_ok and pf_ok:
            valid_options.append((key, loss))
            print("Result: Option is VALID\n")
        else:
            print("Result: Option is INVALID\n")
            
    if not valid_options:
        print("No valid options found based on the constraints.")
        return

    # Find the valid option with the minimum loss
    best_option_key, min_loss = min(valid_options, key=lambda item: item[1])
    
    best_option_values = options[best_option_key]
    P_es_opt = best_option_values['P_ES']
    Q_es_opt = best_option_values['Q_ES']
    
    print("--- Conclusion ---")
    print("The valid options are those that satisfy both E-STATCOM and Power Factor constraints.")
    print(f"Comparing the valid options {valid_options} to find the one with minimum loss.")
    print(f"The minimum loss is {min_loss} MW, which corresponds to option {best_option_key}.")
    print("\nFinal Answer:")
    print("The optimized real and reactive power output of the E-STATCOM are:")
    print(f"P_ES = {P_es_opt} MW")
    print(f"Q_ES = {Q_es_opt} MVAR")
    print("with a total system power loss of:")
    print(f"Loss = {min_loss} MW")

solve_wind_park_optimization()