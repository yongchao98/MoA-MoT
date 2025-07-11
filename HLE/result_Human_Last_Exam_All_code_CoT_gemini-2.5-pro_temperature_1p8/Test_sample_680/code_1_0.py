import math

def solve_tower_placement():
    """
    Formulates and solves the tower placement problem to maximize coverage.
    """
    # Problem Constants
    total_area_km2 = 12 * 11
    total_budget_usd = 45000
    
    # Tower B1 properties
    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * (radius_b1 ** 2)
    # Area of bounding box for packing approximation
    packed_area_b1 = (2 * radius_b1) ** 2 
    
    # Tower B2 properties
    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * (radius_b2 ** 2)
    # Area of bounding box for packing approximation
    packed_area_b2 = (2 * radius_b2) ** 2

    # Variables to store the optimal solution
    best_b1 = 0
    best_b2 = 0
    max_coverage_factor = -1

    # Determine the maximum possible number of B2 towers based on constraints
    max_b2_from_budget = total_budget_usd // cost_b2
    max_b2_from_area = total_area_km2 // packed_area_b2
    max_b2 = min(max_b2_from_budget, max_b2_from_area)

    # Iterate through all possible numbers of b2 towers
    for b2 in range(int(max_b2) + 1):
        # For the current number of b2, find the max possible b1
        remaining_budget = total_budget_usd - b2 * cost_b2
        remaining_area = total_area_km2 - b2 * packed_area_b2
        
        max_b1_from_budget = remaining_budget // cost_b1
        max_b1_from_area = remaining_area // packed_area_b1
        
        # The number of b1 towers is limited by the stricter of the two constraints
        b1 = int(min(max_b1_from_budget, max_b1_from_area))

        # Objective function is proportional to b1 + 4*b2
        current_coverage_factor = b1 + 4 * b2

        # If this combination gives better coverage, store it as the new best
        if current_coverage_factor > max_coverage_factor:
            max_coverage_factor = current_coverage_factor
            best_b1 = b1
            best_b2 = b2

    # Calculate the coverage ratio for the best solution found
    total_coverage_area = best_b1 * area_b1 + best_b2 * area_b2
    coverage_ratio = total_coverage_area / total_area_km2
    coverage_percentage = round(coverage_ratio * 100)

    # Print the result in the specified format b1;b2;c
    # Outputting each number required for the final equation/answer
    print(f"{best_b1};{best_b2};{coverage_percentage}")

solve_tower_placement()
<<<13;5;79>>>