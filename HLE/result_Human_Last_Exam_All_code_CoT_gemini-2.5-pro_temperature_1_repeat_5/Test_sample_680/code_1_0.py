import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem.
    """
    # --- Problem Parameters ---
    RECT_WIDTH = 12  # km
    RECT_HEIGHT = 11 # km
    TOTAL_BUDGET = 45000 # usd

    # Tower B1 specifications
    B1_RADIUS = 1 # km
    B1_COST = 1500 # usd

    # Tower B2 specifications
    B2_RADIUS = 2 # km
    B2_COST = 5000 # usd

    # --- Derived Values ---
    rect_area = RECT_WIDTH * RECT_HEIGHT
    # The objective is to maximize the coverage area, which is proportional to
    # b1 * r1^2 + b2 * r2^2. We call this the 'coverage factor'.
    b1_coverage_factor_unit = B1_RADIUS**2
    b2_coverage_factor_unit = B2_RADIUS**2

    # --- Optimization ---
    # We will iterate through all possible solutions to find the best one.
    best_b1 = 0
    best_b2 = 0
    max_coverage_factor = -1

    # Iterate through all possible numbers of B2 towers within the budget
    max_b2_by_budget = TOTAL_BUDGET // B2_COST
    for b2 in range(max_b2_by_budget + 1):
        # For a given number of b2 towers, find the max b1 towers allowed by budget
        remaining_budget = TOTAL_BUDGET - (b2 * B2_COST)
        b1 = remaining_budget // B1_COST

        # Now we have a candidate solution (b1, b2).

        # Constraint 1: Budget (already handled by the loop structure)

        # Constraint 2: The total area of non-overlapping circles must fit in the rectangle
        current_coverage_factor = b1 * b1_coverage_factor_unit + b2 * b2_coverage_factor_unit
        total_tower_area = math.pi * current_coverage_factor
        
        if total_tower_area <= rect_area:
            # This is a feasible solution. Check if it's better than the best one found so far.
            if current_coverage_factor > max_coverage_factor:
                max_coverage_factor = current_coverage_factor
                best_b1 = b1
                best_b2 = b2

    # --- Final Calculation ---
    # The best combination of towers is (best_b1, best_b2)
    # The maximized coverage area
    max_coverage_area = math.pi * max_coverage_factor
    # The coverage ratio as a percentage
    coverage_ratio = (max_coverage_area / rect_area) * 100
    # Rounded to the nearest integer percentage
    c = round(coverage_ratio)

    # Print the result in the required format: b1;b2;c
    # The numbers in the final equation are best_b1, best_b2, and c.
    print(f"{best_b1};{best_b2};{c}")

solve_tower_placement()
<<<0;9;86>>>