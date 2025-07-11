import math

def solve_tower_placement():
    """
    Solves the WiFi tower placement optimization problem.
    """
    # --- Problem Definition ---
    # City and Budget
    city_area = 12 * 11  # 132 sq km
    total_budget = 45000

    # Tower B1 properties
    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * (radius_b1**2)

    # Tower B2 properties
    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * (radius_b2**2)

    # --- Optimization ---
    # We will find the combination of towers (b1, b2) that maximizes the coverage area
    # while respecting the budget and physical area constraints.

    optimal_b1 = 0
    optimal_b2 = 0
    max_coverage = 0

    # Iterate through all possible numbers of B2 towers.
    # The maximum number of B2 towers is limited by the budget.
    max_b2_possible = total_budget // cost_b2

    for b2 in range(max_b2_possible + 1):
        # For a given number of b2 towers, find the maximum b1 towers we can afford.
        remaining_budget = total_budget - (b2 * cost_b2)
        b1 = remaining_budget // cost_b1

        # Now we have a candidate solution (b1, b2) that meets the budget constraint.
        # We must also check the area constraint (no overlap).
        current_coverage = b1 * area_b1 + b2 * area_b2

        # Check if this solution is valid (fits in the city area) and is the best so far.
        if current_coverage <= city_area:
            if current_coverage > max_coverage:
                max_coverage = current_coverage
                optimal_b1 = b1
                optimal_b2 = b2
    
    # --- Final Result Calculation ---
    # Calculate the coverage ratio as a percentage
    coverage_ratio = (max_coverage / city_area) * 100
    coverage_percent = round(coverage_ratio)

    # Print the values that form the final answer
    print(f"Optimal number of B1 towers (b1): {optimal_b1}")
    print(f"Optimal number of B2 towers (b2): {optimal_b2}")
    print(f"Maximum coverage achieved: {max_coverage:.2f} sq km")
    print(f"Total city area: {city_area} sq km")
    print(f"Coverage Ratio (c): {coverage_percent}%")
    
    # Final output in the required format "b1;b2;c"
    final_answer_string = f"{optimal_b1};{optimal_b2};{coverage_percent}"
    print("\nFormatted Answer:")
    print(final_answer_string)


solve_tower_placement()
<<<0;9;86>>>