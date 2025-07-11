import math

def solve_wifi_optimization():
    """
    Solves the WiFi tower optimization problem to maximize coverage.
    """
    # --- Problem Definition ---
    # City and Budget
    city_area = 12 * 11
    total_budget = 45000

    # Tower B1 specifications
    cost1 = 1500
    radius1 = 1
    area1 = math.pi * radius1**2
    # Objective function weight for b1 is proportional to its area: 1
    weight1 = radius1**2 

    # Tower B2 specifications
    cost2 = 5000
    radius2 = 2
    area2 = math.pi * radius2**2
    # Objective function weight for b2 is proportional to its area: 4
    weight2 = radius2**2

    # --- Find Optimal Solution ---
    # We will iterate through all possible combinations of b1 and b2
    # that are within the constraints and find the one that maximizes coverage.

    best_b1 = 0
    best_b2 = 0
    max_objective_value = -1

    # Determine search range for b2 (from budget: 10*b2 <= 90 => b2 <= 9)
    max_b2 = total_budget // cost2

    # Determine search range for b1 (from budget: 3*b1 <= 90 => b1 <= 30)
    max_b1 = total_budget // cost1
    
    # Iterate through possible numbers of towers
    for b2 in range(max_b2 + 1):
        for b1 in range(max_b1 + 1):
            # Check budget constraint
            current_cost = cost1 * b1 + cost2 * b2
            if current_cost > total_budget:
                continue # This combination is too expensive

            # Check non-overlapping area constraint
            total_weighted_area = weight1 * b1 + weight2 * b2
            if total_weighted_area * math.pi > city_area:
                continue

            # This is a valid combination, check if it's the best so far
            if total_weighted_area > max_objective_value:
                max_objective_value = total_weighted_area
                best_b1 = b1
                best_b2 = b2

    # --- Calculate Final Results ---
    # Total coverage area with the optimal number of towers
    final_coverage_area = best_b1 * area1 + best_b2 * area2
    
    # Coverage ratio as a percentage
    coverage_ratio = (final_coverage_area / city_area) * 100
    coverage_percentage = round(coverage_ratio)
    
    # Print the result in the format b1;b2;c
    # Final equation parts: Number of B1 towers, Number of B2 towers, Coverage Ratio
    print(f"Optimal number of B1 towers: {best_b1}")
    print(f"Optimal number of B2 towers: {best_b2}")
    print(f"Maximized coverage ratio: {coverage_percentage}%")
    
    # Final output string
    final_answer = f"{best_b1};{best_b2};{coverage_percentage}"
    print("\nFormatted Answer:")
    print(final_answer)

solve_wifi_optimization()