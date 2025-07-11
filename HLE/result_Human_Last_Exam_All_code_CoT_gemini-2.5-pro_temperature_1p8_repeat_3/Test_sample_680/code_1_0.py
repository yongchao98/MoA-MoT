import math

def solve_tower_optimization():
    """
    Formulates and solves the tower placement optimization problem.
    """
    # --- Problem Definition ---
    TOTAL_BUDGET = 45000
    CITY_AREA = 12 * 11  # 132 sq km

    # Tower B1 properties
    B1_COST = 1500
    B1_RADIUS = 1
    B1_COVERAGE_AREA = math.pi * B1_RADIUS**2
    B1_BOUNDING_BOX_AREA = (2 * B1_RADIUS)**2  # Area = 4

    # Tower B2 properties
    B2_COST = 5000
    B2_RADIUS = 2
    B2_COVERAGE_AREA = math.pi * B2_RADIUS**2
    B2_BOUNDING_BOX_AREA = (2 * B2_RADIUS)**2  # Area = 16

    # --- ILP Solver ---
    # Objective: Maximize Z = b1 + 4*b2
    # Constraints:
    # 1. 3*b1 + 10*b2 <= 90 (Budget)
    # 2. b1 + 4*b2 <= 33  (Area)

    best_b1 = 0
    best_b2 = 0
    max_objective_value = -1

    # Determine the absolute maximum b2 possible from constraints
    max_b2_from_budget = TOTAL_BUDGET // B2_COST
    max_b2_from_area = CITY_AREA // B2_BOUNDING_BOX_AREA
    max_b2 = min(max_b2_from_budget, max_b2_from_area)

    # Iterate through each possible number of B2 towers
    for b2 in range(max_b2 + 1):
        # For the current number of b2 towers, find the max b1 possible
        remaining_budget = TOTAL_BUDGET - (b2 * B2_COST)
        remaining_area = CITY_AREA - (b2 * B2_BOUNDING_BOX_AREA)
        
        # Ensure budget and area are non-negative
        if remaining_budget < 0 or remaining_area < 0:
            continue
            
        max_b1_from_budget = remaining_budget // B1_COST
        max_b1_from_area = remaining_area // B1_BOUNDING_BOX_AREA
        
        # The number of b1 towers is the most restrictive of the two
        b1 = min(max_b1_from_budget, max_b1_from_area)

        # Calculate the value of the objective function
        current_objective_value = b1 + 4 * b2

        # If this combination gives a better result, store it
        if current_objective_value > max_objective_value:
            max_objective_value = current_objective_value
            best_b1 = b1
            best_b2 = b2

    # --- Calculate Final Coverage ---
    optimal_b1 = best_b1
    optimal_b2 = best_b2

    # Calculate total coverage area with the optimal number of towers
    total_coverage = (optimal_b1 * B1_COVERAGE_AREA) + (optimal_b2 * B2_COVERAGE_AREA)
    
    # Calculate the coverage ratio
    coverage_ratio = total_coverage / CITY_AREA
    
    # Convert to a rounded percentage
    coverage_percentage = round(coverage_ratio * 100)

    # Print the result in the specified format b1;b2;c
    print(f"{optimal_b1};{optimal_b2};{coverage_percentage}")

# Execute the function
solve_tower_optimization()