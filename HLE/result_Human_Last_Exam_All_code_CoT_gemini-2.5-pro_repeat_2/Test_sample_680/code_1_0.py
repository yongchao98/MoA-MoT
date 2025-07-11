import math

def solve_tower_optimization():
    """
    Solves the WiFi tower optimization problem to maximize coverage.
    """
    total_budget = 45000
    rect_area = 12 * 11

    # Using the simplified constraints derived from the plan.
    # 3*b1 + 10*b2 <= 90
    # b1 + 4*b2 <= 38

    best_b1 = -1
    best_b2 = -1
    max_objective_val = -1

    # Iterate through all possible numbers of B2 towers.
    # From 10*b2 <= 90, we know b2 can be at most 9.
    for b2 in range(10):
        # Calculate max b1 based on the budget constraint
        # 3*b1 <= 90 - 10*b2  => b1 <= (90 - 10*b2) / 3
        b1_from_budget = math.floor((90 - 10 * b2) / 3)

        # Calculate max b1 based on the placement constraint
        # b1 <= 38 - 4*b2
        b1_from_packing = 38 - 4 * b2
        
        # b1 must satisfy both constraints, so we take the minimum.
        # Also b1 must be non-negative.
        current_b1 = max(0, min(b1_from_budget, b1_from_packing))

        # Calculate the objective function value
        objective_val = current_b1 + 4 * b2

        # If this combination gives a better result, store it
        if objective_val > max_objective_val:
            max_objective_val = objective_val
            best_b1 = current_b1
            best_b2 = b2

    # Now we have the optimal number of towers, b1 and b2.
    # Let's calculate the final coverage ratio.
    
    # Total coverage area = b1 * pi * r1^2 + b2 * pi * r2^2
    total_coverage_area = (best_b1 * math.pi * 1**2) + (best_b2 * math.pi * 2**2)
    
    # Coverage ratio as a percentage
    coverage_ratio = (total_coverage_area / rect_area) * 100
    
    # Rounded to the nearest integer percentage
    c = round(coverage_ratio)

    # Print the result in the required format
    print(f"{best_b1};{best_b2};{c}")

solve_tower_optimization()