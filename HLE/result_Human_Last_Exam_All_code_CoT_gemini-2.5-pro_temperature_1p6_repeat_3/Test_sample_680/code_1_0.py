import math

def solve_tower_placement():
    """
    Solves the WiFi tower placement optimization problem.
    """
    # --- Problem Parameters ---
    city_area = 12 * 11  # km^2
    total_budget = 45000  # usd
    
    b1_radius = 1  # km
    b1_cost = 1500  # usd
    
    b2_radius = 2  # km
    b2_cost = 5000  # usd

    # --- Initialization for search ---
    best_b1 = 0
    best_b2 = 0
    # The objective is to maximize area, equivalent to maximizing b1 + 4*b2
    max_objective_value = -1

    # --- Iterate through all possible numbers of towers ---
    # We prioritize b2 as it contributes more to the objective function (4*b2 vs 1*b1).
    # The maximum number of B2 towers is limited by the budget.
    max_b2_possible = int(total_budget / b2_cost)
    
    for b2 in range(max_b2_possible, -1, -1):
        # For a given number of b2 towers, find the maximum b1 possible.
        
        # 1. Constraint by budget:
        remaining_budget = total_budget - b2 * b2_cost
        max_b1_from_budget = int(remaining_budget / b1_cost)
        
        # 2. Constraint by area (non-overlapping):
        # b1 * pi*1^2 + b2 * pi*2^2 <= city_area
        # b1 <= city_area/pi - 4*b2
        max_b1_from_area = int((city_area / math.pi) - 4 * b2)

        # The number of b1 towers must satisfy both constraints and be non-negative.
        b1 = max(0, min(max_b1_from_budget, max_b1_from_area))
        
        # Calculate the objective value for this combination
        current_objective_value = b1 + 4 * b2
        
        # If this combination is better, store it as the best so far
        if current_objective_value > max_objective_value:
            max_objective_value = current_objective_value
            best_b1 = b1
            best_b2 = b2

    # --- Calculate Final Results ---
    # Optimal number of towers
    final_b1 = best_b1
    final_b2 = best_b2

    # Area covered by the optimal number of towers
    b1_area = math.pi * (b1_radius ** 2)
    b2_area = math.pi * (b2_radius ** 2)
    total_coverage_area = final_b1 * b1_area + final_b2 * b2_area
    
    # Coverage ratio (as a percentage, rounded to the nearest integer)
    coverage_ratio = (total_coverage_area / city_area) * 100
    final_c = round(coverage_ratio)

    # --- Print the output in the specified format ---
    print(f"{final_b1};{final_b2};{final_c}")

solve_tower_placement()