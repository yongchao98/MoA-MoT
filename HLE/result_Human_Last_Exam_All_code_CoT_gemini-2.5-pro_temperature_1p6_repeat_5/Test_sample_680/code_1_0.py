import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem.
    """
    # --- Problem Parameters ---
    city_width = 12  # km
    city_length = 11  # km
    total_budget = 45000  # usd

    # Tower B1 parameters
    cost_b1 = 1500
    radius_b1 = 1.0
    area_b1 = math.pi * radius_b1**2

    # Tower B2 parameters
    cost_b2 = 5000
    radius_b2 = 2.0
    area_b2 = math.pi * radius_b2**2

    city_area = city_width * city_length

    # --- Optimization ---
    # We will search for the optimal number of towers b1 and b2
    best_b1 = 0
    best_b2 = 0
    max_coverage = 0.0

    # Determine the maximum possible number of each tower based on budget alone
    max_b1_limit = int(total_budget / cost_b1)
    max_b2_limit = int(total_budget / cost_b2)

    # Iterate through all possible combinations of b1 and b2
    for b2 in range(max_b2_limit + 1):
        for b1 in range(max_b1_limit + 1):
            
            current_cost = b1 * cost_b1 + b2 * cost_b2
            current_coverage = b1 * area_b1 + b2 * area_b2
            
            # Check if the current combination meets the constraints
            if current_cost <= total_budget and current_coverage <= city_area:
                # If it's better than what we've found so far, update our best solution
                if current_coverage > max_coverage:
                    max_coverage = current_coverage
                    best_b1 = b1
                    best_b2 = b2
    
    # --- Final Result Calculation ---
    final_cost = best_b1 * cost_b1 + best_b2 * cost_b2
    final_coverage_area = max_coverage
    coverage_ratio = final_coverage_area / city_area
    coverage_percentage = round(coverage_ratio * 100)

    # --- Print the detailed solution ---
    print("Optimization Problem Formulation:")
    print(f"Maximize Coverage = b1 * {area_b1:.2f} + b2 * {area_b2:.2f}")
    print("Subject to:")
    print(f"1. Budget: b1 * {cost_b1} + b2 * {cost_b2} <= {total_budget}")
    print(f"2. Area:   b1 * {area_b1:.2f} + b2 * {area_b2:.2f} <= {city_area}")
    print("\n--- Optimal Solution ---")
    print(f"The optimal number of towers is {best_b1} of type B1 and {best_b2} of type B2.\n")
    
    print("Final Equation Details:")
    # Using the best numbers to show the final equation values
    print(f"Coverage = ({best_b1} * {area_b1:.2f}) + ({best_b2} * {area_b2:.2f}) = {final_coverage_area:.2f} sq km")
    print(f"Cost     = ({best_b1} * {cost_b1}) + ({best_b2} * {cost_b2}) = {final_cost} <= {total_budget} USD")
    print(f"Coverage Ratio = {final_coverage_area:.2f} / {city_area} = {coverage_ratio:.4f} => {coverage_percentage}%\n")

    # --- Final formatted output ---
    final_answer = f"{best_b1};{best_b2};{coverage_percentage}"
    print(f"Final Answer Format: {final_answer}")


solve_tower_placement()
<<<0;9;86>>>