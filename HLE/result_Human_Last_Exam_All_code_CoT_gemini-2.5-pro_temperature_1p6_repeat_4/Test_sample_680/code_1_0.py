import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # --- Problem Definition ---
    city_area_km2 = 12 * 11
    total_budget_usd = 45000

    # Tower B1 properties
    cost_b1 = 1500
    area_b1 = math.pi * (1 ** 2)

    # Tower B2 properties
    cost_b2 = 5000
    area_b2 = math.pi * (2 ** 2)

    # --- Optimization ---
    # We will iterate through all possible numbers of the more expensive tower (B2)
    # and for each, calculate the maximum number of B1 towers we can add.

    best_b1 = 0
    best_b2 = 0
    max_coverage_area = 0.0

    # The maximum number of B2 towers is limited by the total budget.
    max_b2_by_budget = total_budget_usd // cost_b2

    for b2 in range(max_b2_by_budget + 1):
        # Calculate remaining budget and area after placing b2 towers
        remaining_budget = total_budget_usd - (b2 * cost_b2)
        remaining_area = city_area_km2 - (b2 * area_b2)

        # If remaining area is negative, this b2 count is not physically possible.
        if remaining_area < 0:
            continue
            
        # Calculate max B1 towers possible with the remaining resources
        num_b1_by_budget = remaining_budget // cost_b1
        num_b1_by_area = remaining_area // area_b1
        
        # We can only have the minimum of what we can afford and what can fit
        b1 = int(min(num_b1_by_budget, num_b1_by_area))
        
        # Calculate the total coverage for this b1, b2 combination
        current_coverage_area = (b1 * area_b1) + (b2 * area_b2)

        # If this combination is better than the best one found so far, update it.
        if current_coverage_area > max_coverage_area:
            max_coverage_area = current_coverage_area
            best_b1 = b1
            best_b2 = b2

    # --- Final Result ---
    # Calculate the final coverage ratio as a rounded percentage
    coverage_ratio_percent = round((max_coverage_area / city_area_km2) * 100)
    
    # Print the result in the format b1;b2;c
    print(f"{best_b1};{best_b2};{coverage_ratio_percent}")

# Run the solver
solve_tower_placement()
<<<0;9;86>>>