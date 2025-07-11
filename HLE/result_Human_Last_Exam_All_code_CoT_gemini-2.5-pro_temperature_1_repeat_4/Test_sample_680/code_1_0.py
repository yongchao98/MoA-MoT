import math

def solve_tower_optimization():
    """
    Formulates and solves the tower placement problem to maximize coverage.

    This function iterates through all feasible combinations of towers b1 and b2
    to find the one that maximizes the objective function Z = b1 + 4*b2,
    subject to budget and area constraints. It then calculates the coverage
    ratio and prints the final result.
    """
    # City and budget parameters
    city_area = 12 * 11  # km^2
    total_budget = 45000

    # Tower parameters
    cost_b1 = 1500
    cost_b2 = 5000
    radius_b1 = 1
    radius_b2 = 2
    
    # Variables to store the best solution found
    best_b1 = 0
    best_b2 = 0
    max_objective_value = -1

    # Determine the search range for b2. The number of b2 towers cannot exceed
    # the budget or area constraints.
    # Budget limit: 10*b2 <= 90 => b2 <= 9
    # Area limit: 4*b2 <= 42 => b2 <= 10
    # So, the maximum possible b2 is 9.
    max_b2 = total_budget // cost_b2

    # Iterate through all possible numbers of B2 towers
    for b2 in range(max_b2 + 1):
        # For a given b2, find the maximum number of b1 towers possible
        # Budget constraint: 3*b1 <= 90 - 10*b2
        remaining_budget_for_b1 = total_budget - b2 * cost_b2
        max_b1_by_budget = remaining_budget_for_b1 // cost_b1
        
        # Area constraint: b1 <= 42 - 4*b2
        max_b1_by_area = 42 - 4 * b2
        
        # The number of b1 towers is limited by the smaller of the two constraints
        if max_b1_by_area < 0: # Cannot place any more b1 towers
            max_b1 = 0
        else:
            max_b1 = min(max_b1_by_budget, max_b1_by_area)
            
        # We only need to check this maximal b1, as it contributes positively to the objective
        b1 = max_b1
        
        # Check if this combination is valid (it should be by construction, but we double-check)
        budget_ok = (3 * b1 + 10 * b2) <= 90
        area_ok = (b1 + 4 * b2) <= 42

        if budget_ok and area_ok:
            # Calculate the objective function value
            current_objective_value = b1 + 4 * b2
            
            # If this solution is better than the previous best, update it
            if current_objective_value > max_objective_value:
                max_objective_value = current_objective_value
                best_b1 = b1
                best_b2 = b2
    
    # Calculate the final coverage ratio
    # Total coverage area = b1 * pi * r1^2 + b2 * pi * r2^2
    coverage_area_b1 = best_b1 * math.pi * (radius_b1**2)
    coverage_area_b2 = best_b2 * math.pi * (radius_b2**2)
    total_coverage_area = coverage_area_b1 + coverage_area_b2
    
    # Coverage ratio (as a percentage, rounded)
    coverage_ratio = round((total_coverage_area / city_area) * 100)

    # Print the result in the specified format: b1;b2;c
    # The final equation is the result string itself. Here are the numbers that form it:
    print(f"Optimal number of B1 towers (b1): {best_b1}")
    print(f"Optimal number of B2 towers (b2): {best_b2}")
    print(f"Maximized coverage ratio (c): {coverage_ratio}%")
    print("\nFinal Answer:")
    print(f"{best_b1};{best_b2};{coverage_ratio}")


solve_tower_optimization()
<<<0;9;86>>>