import math

def solve_tower_placement():
    """
    Solves the WiFi tower placement optimization problem.
    """
    # --- Problem Parameters ---
    total_budget = 45000
    city_area = 12 * 11

    # B1 Tower parameters
    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * (radius_b1 ** 2)

    # B2 Tower parameters
    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * (radius_b2 ** 2)

    # --- Solving the Optimization Problem ---
    best_b1 = 0
    best_b2 = 0
    max_objective_value = -1

    # The number of B2 towers (b2) is most constrained.
    # From budget: 5000*b2 <= 45000 => b2 <= 9
    # We iterate through all possible numbers of b2 towers.
    max_b2_by_budget = total_budget // cost_b2
    
    for b2 in range(max_b2_by_budget + 1):
        remaining_budget = total_budget - (b2 * cost_b2)
        
        # Calculate the maximum number of b1 towers with the remaining budget
        b1 = remaining_budget // cost_b1

        # Check against the area constraint: b1 + 4*b2 <= 42
        if (b1 + 4 * b2) <= 42:
            # Calculate the objective function value: Z = b1 + 4*b2
            current_objective_value = b1 + 4 * b2
            
            # If this solution is better than the previous best, update it
            if current_objective_value > max_objective_value:
                max_objective_value = current_objective_value
                best_b1 = b1
                best_b2 = b2

    # --- Calculate Final Results ---
    final_b1 = best_b1
    final_b2 = best_b2
    
    # Total covered area for the optimal solution
    total_coverage_area = (final_b1 * area_b1) + (final_b2 * area_b2)
    
    # Coverage ratio as a percentage
    coverage_ratio = (total_coverage_area / city_area) * 100
    coverage_percentage = round(coverage_ratio)

    # --- Print the output as requested ---
    print(f"Optimal Number of B1 towers (b1): {final_b1}")
    print(f"Optimal Number of B2 towers (b2): {final_b2}")
    print(f"Coverage Ratio (c): {coverage_percentage}%")
    
    # Final answer in the specified format b1;b2;c
    final_answer = f"{final_b1};{final_b2};{coverage_percentage}"
    print("\nFormatted Answer:")
    print(final_answer)


solve_tower_placement()
<<<0;9;86>>>