import math

def solve_tower_placement():
    """
    Formulates and solves the tower placement optimization problem.
    """
    # --- 1. Define Problem Parameters ---
    budget_total = 45000
    area_width = 12
    area_length = 11
    total_area = area_width * area_length

    # B1 tower specs
    b1_radius = 1
    b1_cost = 1500
    b1_coverage_area = math.pi * b1_radius**2
    b1_pack_area = (2 * b1_radius)**2  # Effective area in a grid pack (2x2)

    # B2 tower specs
    b2_radius = 2
    b2_cost = 5000
    b2_coverage_area = math.pi * b2_radius**2
    b2_pack_area = (2 * b2_radius)**2  # Effective area in a grid pack (4x4)

    # --- 2. Solve the Integer Program ---
    best_b1 = -1
    best_b2 = -1
    max_coverage_objective = -1
    min_cost_at_max = float('inf')

    # Determine the search range for b2, the more impactful variable
    max_b2_from_budget = math.floor(budget_total / b2_cost)
    max_b2_from_packing = math.floor(total_area / b2_pack_area)
    max_b2 = min(max_b2_from_budget, max_b2_from_packing)

    # Iterate through all valid numbers of towers
    for num_b2 in range(max_b2 + 1):
        # Given num_b2, find the max possible num_b1 allowed by constraints
        remaining_budget = budget_total - num_b2 * b2_cost
        remaining_pack_area = total_area - num_b2 * b2_pack_area

        if remaining_budget < 0 or remaining_pack_area < 0:
            continue

        max_b1_from_budget = math.floor(remaining_budget / b1_cost)
        max_b1_from_packing = math.floor(remaining_pack_area / b1_pack_area)
        
        num_b1 = min(max_b1_from_budget, max_b1_from_packing)
        
        # Calculate the value of the objective function (proportional to coverage)
        current_coverage_objective = num_b1 + 4 * num_b2
        
        # Check if this solution is better than the current best
        if current_coverage_objective > max_coverage_objective:
            max_coverage_objective = current_coverage_objective
            best_b1 = num_b1
            best_b2 = num_b2
            min_cost_at_max = num_b1 * b1_cost + num_b2 * b2_cost
        elif current_coverage_objective == max_coverage_objective:
            # Tie-breaker: choose the solution with lower cost
            current_cost = num_b1 * b1_cost + num_b2 * b2_cost
            if current_cost < min_cost_at_max:
                best_b1 = num_b1
                best_b2 = num_b2
                min_cost_at_max = current_cost

    # --- 3. Calculate Final Results ---
    b1_final = best_b1
    b2_final = best_b2

    # Final calculations using the optimal numbers
    total_covered_area = b1_final * b1_coverage_area + b2_final * b2_coverage_area
    coverage_ratio = total_covered_area / total_area
    coverage_percentage = round(coverage_ratio * 100)
    
    # --- 4. Print detailed output ---
    print(f"The optimization goal is to maximize the expression 'b1 + 4*b2' subject to budget and packing constraints.")
    print(f"\nOptimal solution found: {b1_final} B1 towers and {b2_final} B2 towers.\n")
    
    final_cost = b1_final * b1_cost + b2_final * b2_cost
    final_pack = b1_final * b1_pack_area + b2_final * b2_pack_area

    print("Verification of the final solution against constraints:")
    print(f"Budget Equation: {b1_final} * {b1_cost} + {b2_final} * {b2_cost} = {final_cost} (<= {budget_total})")
    print(f"Packing Equation: {b1_final} * {b1_pack_area} + {b2_final} * {b2_pack_area} = {final_pack} (<= {total_area})\n")

    print(f"This configuration results in a total coverage of {total_covered_area:.2f} sq km.")
    print(f"The coverage ratio for the {total_area} sq km city is {coverage_ratio:.2%}.")
    print(f"The final rounded coverage percentage (c) is {coverage_percentage}%.")

    # Final formatted answer
    final_answer_string = f"{b1_final};{b2_final};{coverage_percentage}"
    print(f"\n<<<ANSER>>>\n{final_answer_string}")

solve_tower_placement()