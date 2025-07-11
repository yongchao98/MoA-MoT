def calculate_t_gate_overhead():
    """
    Calculates the approximate number of non-Clifford gates required for
    two scenarios in topological quantum computing based on surface codes.
    """
    # --- Input Parameters ---
    p_phys = 0.01  # Physical gate error rate (1%)
    d1 = 3         # Distance for Scenario 1 ("simulation")
    d2 = 5         # Distance for Scenario 2 ("implementation")

    # Magic state distillation protocol parameters (15-to-1)
    distill_input_count = 15
    # Error reduction formula: p_out = C * p_in^3. For the 15-to-1 protocol, C is approximately 35.
    distill_error_coeff = 35
    distill_error_power = 3

    print("This script calculates the number of non-Clifford T-gates needed to produce one high-fidelity logical T-gate.")
    print(f"Assumption: A physical gate error rate of {p_phys*100}%.")
    print("Assumption: The 15-to-1 magic state distillation protocol is used.")
    print("-" * 70)

    # --- Scenario 1: distance-3 code ---
    print(f"Scenario 1: Simulation on a 2D Surface Code (distance d={d1})")
    
    # Step 1: Calculate the target logical error rate for Clifford gates
    p_target_1_exp = (d1 + 1) / 2
    p_target_1 = p_phys ** p_target_1_exp
    print(f"Step 1: Calculate the target error rate for a d={d1} code.")
    print(f"   - Target Logical Error (p_L) ≈ p_phys^((d+1)/2) = {p_phys}^({int(p_target_1_exp)}) = {p_target_1:.2e}")

    # Step 2: Calculate the error after one level of distillation
    levels_1 = 1
    p_distilled_1 = distill_error_coeff * (p_phys ** distill_error_power)
    print("Step 2: Calculate the T-gate error after one level of distillation.")
    print(f"   - p_out = {distill_error_coeff} * p_in^{distill_error_power} = {distill_error_coeff} * {p_phys}^{distill_error_power} = {p_distilled_1:.2e}")
    
    # Step 3: Compare and conclude
    print(f"Step 3: Compare distilled error ({p_distilled_1:.2e}) with target ({p_target_1:.2e}).")
    if p_distilled_1 <= p_target_1:
        print("   - Conclusion: One level of distillation is sufficient as the distilled error is lower than the target.")
    else:
        # This case is not reached with the given numbers
        print("   - Conclusion: More distillation levels are needed.")
        levels_1 = 0 # Mark as failed for safety, though it won't happen here.

    num_gates_1 = distill_input_count ** levels_1
    print("Step 4: Calculate the total physical non-Clifford gates required.")
    print(f"   - Equation for d={d1}: {distill_input_count} ^ {levels_1}")
    print(f"   - Result for d={d1}: {num_gates_1} gates")
    
    print("-" * 70)

    # --- Scenario 2: distance-5 code ---
    print(f"Scenario 2: Implementation on a 2D Surface Code (distance d={d2})")
    
    # Step 1: Calculate target logical error rate
    p_target_2_exp = (d2 + 1) / 2
    p_target_2 = p_phys ** p_target_2_exp
    print(f"Step 1: Calculate the target error rate for a d={d2} code.")
    print(f"   - Target Logical Error (p_L) ≈ p_phys^((d+1)/2) = {p_phys}^({int(p_target_2_exp)}) = {p_target_2:.2e}")
    
    # Step 2: Check error after one level of distillation
    print(f"Step 2: Check if one level of distillation is sufficient.")
    print(f"   - Error after one level is {p_distilled_1:.2e}. Target is {p_target_2:.2e}.")
    if p_distilled_1 <= p_target_2:
        print(f"   - Conclusion: One level of distillation is sufficient.")
        levels_2 = 1
    else:
        print(f"   - Conclusion: One level is not enough. A second level is needed.")
        levels_2 = 2
        
        # Step 3: Calculate error after a second level of distillation
        p_distilled_2 = distill_error_coeff * (p_distilled_1 ** distill_error_power)
        print("Step 3: Calculate the T-gate error after a second level of distillation.")
        print(f"   - p_out_L2 = {distill_error_coeff} * p_out_L1^{distill_error_power} = {distill_error_coeff} * ({p_distilled_1:.2e})^{distill_error_power} = {p_distilled_2:.2e}")
        if p_distilled_2 <= p_target_2:
             print(f"   - The new error ({p_distilled_2:.2e}) is now below the target ({p_target_2:.2e}). Two levels are sufficient.")
        else:
             print("   - More levels would be needed.")
             levels_2 = 0 # Mark as failed for safety

    num_gates_2 = distill_input_count ** levels_2
    print(f"Step 4: Calculate the total physical non-Clifford gates required.")
    if levels_2 == 1:
        print(f"   - Equation for d={d2}: {distill_input_count} ^ {levels_2}")
    else:
        print(f"   - Equation for d={d2}: {distill_input_count} ^ {levels_2} = {distill_input_count} * {distill_input_count}")
    print(f"   - Result for d={d2}: {num_gates_2} gates")
    
    print("-" * 70)
    
    # --- Final Summary ---
    total_gates = num_gates_1 + num_gates_2
    print("Final Calculation: Total non-Clifford gates for both tasks sequentially.")
    print(f"Total Gates = (Gates for Scenario 1) + (Gates for Scenario 2)")
    print(f"Final Equation: {num_gates_1} + {num_gates_2} = {total_gates}")


if __name__ == '__main__':
    calculate_t_gate_overhead()