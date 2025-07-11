import math

def calculate_gate_costs():
    """
    Calculates the approximate number of non-Clifford T-gates required for
    fault-tolerant operations on distance-3 and distance-5 surface codes.
    """
    # 1. Define initial parameters
    p_phys = 0.01  # Physical gate error rate
    distillation_input_size = 15 # For a 15-to-1 distillation protocol

    print(f"Initial physical gate error rate (p_phys): {p_phys}")
    print("Using the 15-to-1 magic state distillation protocol.\n")

    # --- Part 1: Cost for Distance-3 Code ---
    d1 = 3
    print(f"--- Calculating cost for distance-3 code (d={d1}) ---")

    # 2. Calculate the target logical error rate for d=3
    # Formula: p_L ≈ 10 * p_phys^((d+1)/2)
    p_logical_d3 = 10 * (p_phys**((d1 + 1) / 2))
    print(f"Target logical error rate (p_L) for d={d1} is approximately: {p_logical_d3:.2e}")

    # 3. Calculate error after one round of distillation
    # Formula: p_out ≈ 35 * p_in^3
    p_distilled_r1 = 35 * (p_phys**3)
    print(f"Error rate after 1 round of distillation: {p_distilled_r1:.2e}")

    # 4. Determine cost for d=3
    if p_distilled_r1 < p_logical_d3:
        cost_d3 = distillation_input_size
        print(f"One round is sufficient for d={d1}. Cost = {cost_d3} T-gates.")
    else:
        # This case is not expected for d=3 with these parameters
        cost_d3 = distillation_input_size * distillation_input_size
        print(f"One round is NOT sufficient for d={d1}. Needing two rounds. Cost = {cost_d3} T-gates.")

    print("\n" + "="*50 + "\n")

    # --- Part 2: Cost for Distance-5 Code ---
    d2 = 5
    print(f"--- Calculating cost for distance-5 code (d={d2}) ---")

    # 5. Calculate the target logical error rate for d=5
    p_logical_d5 = 10 * (p_phys**((d2 + 1) / 2))
    print(f"Target logical error rate (p_L) for d={d2} is approximately: {p_logical_d5:.2e}")
    print(f"Error rate after 1 round of distillation (from above): {p_distilled_r1:.2e}")

    # 6. Determine cost for d=5
    if p_distilled_r1 < p_logical_d5:
        cost_d5 = distillation_input_size
        print(f"One round is sufficient for d={d2}. Cost = {cost_d5} T-gates.")
    else:
        print(f"One round is NOT sufficient. A second round of distillation is needed.")
        # Error after a second round of distillation
        p_distilled_r2 = 35 * (p_distilled_r1**3)
        cost_d5 = distillation_input_size * distillation_input_size
        print(f"Error rate after 2 rounds of distillation: {p_distilled_r2:.2e}")
        print(f"Cost for two rounds = {distillation_input_size} * {distillation_input_size} = {cost_d5} T-gates.")

    print("\n" + "="*50 + "\n")

    # --- Part 3: Final Summation ---
    total_cost = cost_d3 + cost_d5
    print("--- Total Cost Calculation ---")
    print(f"The total approximate number of non-Clifford gates is the sum of the costs for both implementations.")
    print(f"Total Gates = (Cost for d=3) + (Cost for d=5)")
    print(f"Total Gates = {cost_d3} + {cost_d5} = {total_cost}")

    return total_cost

# Execute the calculation and print the final answer in the required format
final_answer = calculate_gate_costs()
print(f"\n<<<__{final_answer}__>>>")
