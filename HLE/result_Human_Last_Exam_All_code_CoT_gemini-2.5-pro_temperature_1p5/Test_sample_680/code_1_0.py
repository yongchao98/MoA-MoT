import math

def solve_tower_optimization():
    """
    Solves the tower placement optimization problem by iterating through all
    feasible solutions to find the one that maximizes coverage.
    """
    # --- Problem Parameters ---
    total_budget = 45000
    rect_width = 12
    rect_length = 11
    total_area = rect_width * rect_length

    b1_radius = 1
    b1_cost = 1500
    b1_area = math.pi * b1_radius**2
    b1_objective_coeff = 1 # from (b1 * 1^2)

    b2_radius = 2
    b2_cost = 5000
    b2_area = math.pi * b2_radius**2
    b2_objective_coeff = 4 # from (b2 * 2^2)
    
    # --- Optimization ---
    best_b1 = 0
    best_b2 = 0
    max_objective_val = -1

    # Iterate through all possible numbers of b2 towers
    # Max b2 is limited by budget (45000/5000 = 9)
    max_b2_by_budget = total_budget // b2_cost
    
    for b2 in range(max_b2_by_budget + 1):
        # For a given b2, find the maximum possible b1
        
        # 1. Budget constraint on b1
        remaining_budget = total_budget - b2 * b2_cost
        max_b1_by_budget = remaining_budget // b1_cost
        
        # 2. Area constraint on b1
        # b1 * pi + b2 * 4pi <= total_area  => b1 <= total_area/pi - 4*b2
        # Since b1 must be an integer, we take the floor.
        max_b1_by_area = math.floor(total_area / math.pi - b2 * b2_objective_coeff)
        
        # The number of b1 towers must satisfy both constraints
        if max_b1_by_budget < 0 or max_b1_by_area < 0:
            continue
            
        b1 = min(max_b1_by_budget, max_b1_by_area)
        
        # Calculate the objective function value for this combination
        current_objective_val = b1 * b1_objective_coeff + b2 * b2_objective_coeff
        
        # If this solution is better than the best one found so far, update
        if current_objective_val > max_objective_val:
            max_objective_val = current_objective_val
            best_b1 = b1
            best_b2 = b2

    # --- Calculate Final Results ---
    final_b1 = best_b1
    final_b2 = best_b2

    # Calculate the total cost and covered area for the final solution
    final_cost = final_b1 * b1_cost + final_b2 * b2_cost
    final_covered_area = final_b1 * b1_area + final_b2 * b2_area

    # Calculate coverage ratio as a percentage
    coverage_ratio = (final_covered_area / total_area) * 100
    coverage_percentage = round(coverage_ratio)

    # --- Print the final output in the required format ---
    # The final equation components: b1, b2, and the resulting coverage percentage 'c'
    print(f"{final_b1};{final_b2};{coverage_percentage}")

solve_tower_optimization()