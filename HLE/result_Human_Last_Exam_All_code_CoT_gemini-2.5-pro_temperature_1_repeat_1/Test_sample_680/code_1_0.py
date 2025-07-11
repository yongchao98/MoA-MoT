import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem.
    """
    # --- Problem Parameters ---
    city_width = 12
    city_height = 11
    total_budget = 45000

    # B1 Tower
    radius1 = 1
    cost1 = 1500

    # B2 Tower
    radius2 = 2
    cost2 = 5000

    # --- Calculations ---
    city_area = city_width * city_height
    
    best_b1 = 0
    best_b2 = 0
    max_objective_value = -1

    # Iterate through all possible numbers of B2 towers
    # The maximum number of B2 towers is limited by the budget
    max_b2_possible = total_budget // cost2

    for b2 in range(max_b2_possible + 1):
        # For a given number of b2 towers, find the max b1 towers
        remaining_budget = total_budget - (b2 * cost2)
        b1 = remaining_budget // cost1

        # Check the area constraint
        # The total area of circles must not exceed the city area
        total_circle_area = math.pi * (b1 * radius1**2 + b2 * radius2**2)
        
        if total_circle_area <= city_area:
            # Calculate the objective function value (proportional to coverage)
            current_objective_value = b1 * radius1**2 + b2 * radius2**2
            
            # If this is the best solution so far, save it
            if current_objective_value > max_objective_value:
                max_objective_value = current_objective_value
                best_b1 = b1
                best_b2 = b2

    # --- Final Result ---
    # Calculate the total cost and coverage for the best solution
    final_cost = best_b1 * cost1 + best_b2 * cost2
    final_coverage_area = math.pi * max_objective_value
    coverage_ratio = (final_coverage_area / city_area) * 100
    coverage_percentage = round(coverage_ratio)

    # Print the details of the final solution as requested
    print(f"Optimal Solution Found:")
    print(f"Number of B1 towers (b1): {best_b1}")
    print(f"Number of B2 towers (b2): {best_b2}")
    
    print("\nVerification of the final equation values:")
    print(f"Objective Value (b1 + 4*b2): {best_b1} + 4*{best_b2} = {max_objective_value}")
    print(f"Budget Used: {cost1}*{best_b1} + {cost2}*{best_b2} = ${final_cost} (<= ${total_budget})")
    print(f"Area Used: pi*({best_b1}*1^2 + {best_b2}*2^2) = {final_coverage_area:.2f} sq.km (<= {city_area} sq.km)")

    # Print the final answer in the required format
    print("\nFinal Answer (b1;b2;coverage_percentage):")
    print(f"{best_b1};{best_b2};{coverage_percentage}")

solve_tower_placement()