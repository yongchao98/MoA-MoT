import math

def solve_quantum_computing_cost():
    """
    Calculates the approximate number of non-Clifford gates required for universal
    quantum computation on surface codes of distance 3 and 5.
    """
    # Physical gate error rate
    p_physical = 0.01

    # --- Magic State Distillation Protocol Parameters (15-to-1 protocol) ---
    # Number of input states consumed to produce one output state
    distillation_overhead = 15
    # Factor in the error reduction formula: p_out = factor * p_in^3
    distillation_error_factor = 35

    print("This script estimates the number of non-Clifford gates needed for fault-tolerant computation.")
    print(f"Assuming a physical gate error rate of {p_physical*100}% and a 15-to-1 magic state distillation protocol.\n")

    # --- Scenario 1: 2D Surface Code with distance d=3 ---
    d1 = 3
    print(f"--- Scenario 1: Distance-3 Code (d={d1}) ---")

    # Step 1: Calculate the target logical error rate for d=3
    # Formula: p_L ≈ 10 * p^((d+1)/2)
    p_target_1 = 10 * (p_physical**((d1 + 1) / 2))
    print(f"The target logical error rate is approximately {p_target_1:.1e}.")

    # Step 2: Calculate the number of distillation rounds and the gate cost
    cost1 = 1
    num_rounds1 = 0
    current_error1 = p_physical
    while current_error1 > p_target_1:
        num_rounds1 += 1
        cost1 *= distillation_overhead
        p_out = distillation_error_factor * (current_error1**3)
        print(f"After round {num_rounds1} (total cost: {cost1} gates), the error is reduced from {current_error1:.1e} to {p_out:.1e}.")
        current_error1 = p_out
    
    print(f"The final error {current_error1:.1e} is below the target. Required gates for d=3: {cost1}\n")

    # --- Scenario 2: 2D Surface Code with distance d=5 ---
    d2 = 5
    print(f"--- Scenario 2: Distance-5 Code (d={d2}) ---")

    # Step 1: Calculate the target logical error rate for d=5
    p_target_2 = 10 * (p_physical**((d2 + 1) / 2))
    print(f"The target logical error rate is approximately {p_target_2:.1e}.")

    # Step 2: Calculate the number of distillation rounds and the gate cost
    cost2 = 1
    num_rounds2 = 0
    current_error2 = p_physical
    while current_error2 > p_target_2:
        num_rounds2 += 1
        cost2 *= distillation_overhead
        p_out = distillation_error_factor * (current_error2**3)
        print(f"After round {num_rounds2} (total cost: {cost2} gates), the error is reduced from {current_error2:.1e} to {p_out:.1e}.")
        current_error2 = p_out

    print(f"The final error {current_error2:.1e} is below the target. Required gates for d=5: {cost2}\n")

    # --- Final Calculation ---
    total_cost = cost1 + cost2
    print("--- Total Cost ---")
    print(f"Total approximate number of non-Clifford gates = {cost1} + {cost2} = {total_cost}")
    
    return total_cost

# Run the calculation and store the final answer
final_answer = solve_quantum_computing_cost()
# The final answer is wrapped as requested by the user.
# print(f"<<<{final_answer}>>>")