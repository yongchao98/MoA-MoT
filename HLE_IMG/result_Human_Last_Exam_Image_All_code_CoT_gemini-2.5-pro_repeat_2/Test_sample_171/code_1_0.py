import math

def solve_optimization_problem():
    """
    Evaluates the given choices for the E-STATCOM operation against the system constraints.
    """
    # System Parameters
    P_wp = 8.0  # MW
    P_load = 6.0  # MW
    Q_load = 2.0  # MVAR
    S_base = 10.0  # MVA
    V_grid_pu = 1.0 # per unit

    # Line Impedance (assuming given values are on 100 MVA base, converting to 10 MVA base)
    R_line_pu_100MVA = 0.05
    X_line_pu_100MVA = 0.2
    S_base_old = 100.0
    S_base_new = 10.0
    R_line_pu = R_line_pu_100MVA * (S_base_new / S_base_old)
    X_line_pu = X_line_pu_100MVA * (S_base_new / S_base_old)

    # Constraints
    P_ES_max = 4.0  # MW
    Q_ES_max = 3.0  # MVAR
    PF_min = 0.98
    V_dev_max = 0.015  # 1.5%

    # Answer Choices [P_ES (MW), Q_ES (MVAR), Total Loss (MW)]
    choices = {
        'A': [1.5, 1.8, 0.5],
        'B': [3.2, 2.1, 0.4],
        'C': [3.0, 2.0, 0.45],
        'D': [2.5, 2.8, 0.35],
        'E': [3.5, 2.5, 0.5]
    }

    print("Evaluating Answer Choices:\n")
    
    feasible_choices = {}

    for name, values in choices.items():
        P_ES, Q_ES, loss = values
        print(f"--- Checking Choice {name} ---")
        print(f"P_ES = {P_ES} MW, Q_ES = {Q_ES} MVAR, Stated Loss = {loss} MW")

        # --- Constraint Checks ---
        valid = True

        # 1. E-STATCOM Limits
        if not (abs(P_ES) <= P_ES_max and abs(Q_ES) <= Q_ES_max):
            valid = False
            print("Status: FAILED - E-STATCOM power limits exceeded.")
            continue
        
        # Calculate PCC power injection
        P_pcc = P_wp + P_ES
        Q_pcc = Q_ES

        # 2. Power Factor Constraint
        S_pcc = math.sqrt(P_pcc**2 + Q_pcc**2)
        if S_pcc == 0: # Avoid division by zero
            pf = 1.0
        else:
            pf = P_pcc / S_pcc
        
        if pf < PF_min:
            valid = False
            print(f"Power Factor Check: Calculated PF = {pf:.4f}, Required >= {PF_min}")
            print("Status: FAILED - Power factor requirement not met.")
            continue
        else:
            print(f"Power Factor Check: Calculated PF = {pf:.4f}, Required >= {PF_min} -> OK")

        # 3. Voltage Constraint
        P_pcc_pu = P_pcc / S_base
        Q_pcc_pu = Q_pcc / S_base
        
        # V_pcc_pu = V_grid_pu + delta_V_pu
        delta_V_pu = (P_pcc_pu * R_line_pu) + (Q_pcc_pu * X_line_pu)
        V_pcc_pu = V_grid_pu + delta_V_pu
        
        voltage_deviation = abs(V_pcc_pu - V_grid_pu)
        
        if voltage_deviation > V_dev_max:
            valid = False
            print(f"Voltage Check: Calculated V_pcc = {V_pcc_pu:.4f} p.u., Deviation = {voltage_deviation*100:.2f}%")
            print(f"Status: FAILED - Voltage deviation exceeds {V_dev_max*100}%.")
            continue
        else:
            print(f"Voltage Check: Calculated V_pcc = {V_pcc_pu:.4f} p.u., Deviation = {voltage_deviation*100:.2f}% -> OK")

        if valid:
            print("Status: PASSED - All constraints are satisfied.")
            feasible_choices[name] = values
        print("-" * 25)

    print("\n--- Conclusion ---")
    if not feasible_choices:
        print("No feasible choices found based on the problem statement and assumptions.")
    else:
        print("Feasible choices are:", list(feasible_choices.keys()))
        # Find the feasible choice with the minimum stated loss
        best_choice_name = min(feasible_choices, key=lambda k: feasible_choices[k][2])
        best_choice_values = feasible_choices[best_choice_name]
        
        P_ES_opt = best_choice_values[0]
        Q_ES_opt = best_choice_values[1]
        loss_opt = best_choice_values[2]

        print(f"\nThe optimal solution is the feasible choice with the minimum power loss.")
        print(f"The minimum loss among feasible options is {loss_opt} MW, which corresponds to Choice {best_choice_name}.")
        print(f"Final Answer: The optimized real and reactive power output of the E-STATCOM are P_ES = {P_ES_opt} MW and Q_ES = {Q_ES_opt} MVAR, with a total system power loss of {loss_opt} MW.")

solve_optimization_problem()