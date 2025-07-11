import math

def solve_tower_optimization():
    """
    Solves the tower placement optimization problem to maximize coverage.

    This function formulates the problem as an integer optimization task (specifically,
    a knapsack-type problem) to find the ideal number of B1 and B2 towers that
    maximizes area coverage within a given budget.

    The geometric constraint of non-overlapping placement is simplified by assuming
    a valid placement is possible for any chosen combination of towers. The code
    iterates through all feasible combinations to find the optimal solution.
    """

    # --- Problem Definition ---
    city_area = 12 * 11  # km^2
    total_budget = 45000  # USD

    # Tower B1 properties
    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * radius_b1**2
    # The 'value' in our optimization is proportional to the area.
    # We want to maximize n1*area_b1 + n2*area_b2 = pi * (n1*1 + n2*4)
    # So we maximize the value function V = n1 + 4*n2
    value_b1 = 1 

    # Tower B2 properties
    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * radius_b2**2
    value_b2 = 4

    # --- Optimization ---
    best_n1 = 0
    best_n2 = 0
    max_value = 0

    # Iterate through all possible numbers of B2 towers
    max_possible_n2 = total_budget // cost_b2
    for n2 in range(max_possible_n2 + 1):
        # Calculate the remaining budget after placing n2 towers
        remaining_budget = total_budget - (n2 * cost_b2)
        
        # Calculate the maximum number of B1 towers we can afford with the rest
        n1 = remaining_budget // cost_b1
        
        # Calculate the objective function value for this combination
        current_value = (n1 * value_b1) + (n2 * value_b2)
        
        # If this combination gives a better value, store it as the best so far
        if current_value > max_value:
            max_value = current_value
            best_n1 = n1
            best_n2 = n2

    # --- Calculate Final Results ---
    # The optimal number of towers
    b1_final_count = best_n1
    b2_final_count = best_n2

    # Calculate the total area covered by the optimal number of towers
    total_coverage = (b1_final_count * area_b1) + (b2_final_count * area_b2)
    
    # Calculate the coverage ratio and round it to the nearest percentage
    coverage_ratio = total_coverage / city_area
    coverage_percentage = round(coverage_ratio * 100)

    # --- Output the result in the specified format ---
    # The final equation's components are the number of B1 towers, B2 towers, and the coverage percentage
    print(f"{b1_final_count};{b2_final_count};{coverage_percentage}")

solve_tower_optimization()