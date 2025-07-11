import math

def solve_tower_placement():
    """
    Finds the optimal number of B1 and B2 towers to maximize coverage
    within budget and area constraints.
    """
    # --- Problem Parameters ---
    CITY_AREA = 12 * 11  # km^2
    TOTAL_BUDGET = 45000  # usd

    B1_RADIUS = 1  # km
    B1_COST = 1500  # usd
    B1_AREA = math.pi * (B1_RADIUS ** 2)

    B2_RADIUS = 2  # km
    B2_COST = 5000  # usd
    B2_AREA = math.pi * (B2_RADIUS ** 2)

    # --- Optimization ---
    best_b1 = 0
    best_b2 = 0
    max_coverage = 0.0

    # Iterate through all possible numbers of B2 towers
    # The max number of B2 towers is limited by the budget
    max_num_b2 = TOTAL_BUDGET // B2_COST
    for b2 in range(max_num_b2 + 1):
        
        # For each number of B2 towers, calculate the remaining budget for B1 towers
        remaining_budget = TOTAL_BUDGET - (b2 * B2_COST)
        max_num_b1 = remaining_budget // B1_COST
        
        # We only need to check the case with the max number of B1 towers,
        # as that will always yield the highest coverage for a fixed number of B2 towers.
        b1 = max_num_b1

        current_coverage = (b1 * B1_AREA) + (b2 * B2_AREA)

        # Check if this combination is valid (fits within the city area)
        # and if it's better than the best one found so far.
        if current_coverage <= CITY_AREA and current_coverage > max_coverage:
            max_coverage = current_coverage
            best_b1 = b1
            best_b2 = b2
    
    # --- Final Calculation ---
    coverage_ratio = max_coverage / CITY_AREA
    # Round to the nearest percentage
    coverage_percentage = int(round(coverage_ratio * 100))

    # Output the result in the specified format
    print(f"{best_b1};{best_b2};{coverage_percentage}")

solve_tower_placement()