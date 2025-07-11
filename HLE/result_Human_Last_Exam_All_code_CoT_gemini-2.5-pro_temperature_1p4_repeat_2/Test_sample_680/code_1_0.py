import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # Define problem constants
    city_width = 12  # km
    city_length = 11 # km
    total_budget = 45000 # usd

    # Tower B1 properties
    radius_b1 = 1 # km
    cost_b1 = 1500 # usd
    area_b1 = math.pi * radius_b1**2

    # Tower B2 properties
    radius_b2 = 2 # km
    cost_b2 = 5000 # usd
    area_b2 = math.pi * radius_b2**2

    # Total city area
    city_area = city_width * city_length

    # --- Step 1: Find the optimal number of towers ---
    
    # We will iterate through all possible combinations that respect the budget
    # and find the one that maximizes the coverage area.

    best_b1 = 0
    best_b2 = 0
    max_coverage_area = 0

    # Determine the maximum possible number of each tower type based on budget
    max_num_b2 = total_budget // cost_b2
    
    # Iterate through all possible numbers of B2 towers (from max to 0)
    for num_b2 in range(max_num_b2 + 1):
        # Calculate the remaining budget after placing B2 towers
        remaining_budget = total_budget - (num_b2 * cost_b2)
        
        # Calculate the maximum number of B1 towers that can be afforded
        num_b1 = remaining_budget // cost_b1

        # Calculate the total cost and coverage for this combination
        current_cost = num_b1 * cost_b1 + num_b2 * cost_b2
        current_coverage_area = num_b1 * area_b1 + num_b2 * area_b2

        # Check against constraints
        # 1. Budget constraint is implicitly handled by the loops
        # 2. Area constraint (total tower area must be less than city area)
        if current_coverage_area <= city_area:
            # Check if this combination gives better coverage
            if current_coverage_area > max_coverage_area:
                max_coverage_area = current_coverage_area
                best_b1 = num_b1
                best_b2 = num_b2

    # --- Step 2: Calculate and print the results ---
    
    final_cost = best_b1 * cost_b1 + best_b2 * cost_b2
    coverage_ratio_percentage = round((max_coverage_area / city_area) * 100)

    print(f"Optimization Problem Formulation:")
    print(f"Maximize: b1 * {area_b1:.2f} + b2 * {area_2:.2f}")
    print(f"Subject to:")
    print(f"  b1 * {cost_b1} + b2 * {cost_b2} <= {total_budget} (Budget)")
    print(f"  b1 * {area_b1:.2f} + b2 * {area_2:.2f} <= {city_area} (Area)\n")

    print(f"Optimal Solution Found:")
    print(f"- Number of B1 towers (b1): {best_b1}")
    print(f"- Number of B2 towers (b2): {best_b2}\n")

    print(f"Verification:")
    print(f"- Total Cost: {best_b1} * {cost_b1} + {best_b2} * {cost_b2} = ${final_cost:.2f} (<= ${total_budget:.2f})")
    print(f"The Final Equation for Coverage Area:")
    print(f"Maximized Coverage Area: {best_b1} * {area_b1:.2f} km^2 + {best_b2} * {area_2:.2f} km^2 = {max_coverage_area:.2f} km^2\n")
    
    print(f"Coverage Ratio:")
    print(f"Coverage Ratio = ({max_coverage_area:.2f} km^2 / {city_area} km^2) * 100% = {coverage_ratio_percentage}%\n")

    # Final formatted output
    final_answer_string = f"{best_b1};{best_b2};{coverage_ratio_percentage}"
    print("Final Answer String:")
    print(final_answer_string)
    
    # Return final answer in the required format
    print(f"<<<{final_answer_string}>>>")


solve_tower_placement()