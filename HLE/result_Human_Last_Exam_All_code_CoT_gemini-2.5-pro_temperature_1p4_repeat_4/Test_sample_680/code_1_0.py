import math

def solve_wifi_optimization():
    """
    Formulates and solves the WiFi tower placement optimization problem.
    """
    # 1. Define problem parameters
    city_width = 12
    city_length = 11
    city_area = city_width * city_length

    b1_radius = 1
    b1_cost = 1500
    b1_area = math.pi * b1_radius**2
    b1_bounding_square = (2 * b1_radius)**2

    b2_radius = 2
    b2_cost = 5000
    b2_area = math.pi * b2_radius**2
    b2_bounding_square = (2 * b2_radius)**2

    total_budget = 45000

    # 2. Iterate to find the optimal solution
    best_b1 = 0
    best_b2 = 0
    max_objective_value = -1
    highest_cost_for_max_objective = -1

    # Determine search range for b2 from constraints to be efficient
    # From budget: 5000*b2 <= 45000 -> b2 <= 9
    # From packing: 16*b2 <= 132 -> b2 <= 8.25
    max_b2 = int(min(total_budget / b2_cost, city_area / b2_bounding_square))

    for b2 in range(max_b2 + 1):
        # Determine max b1 for the current b2
        limit_from_budget = (total_budget - b2 * b2_cost) / b1_cost
        limit_from_packing = (city_area - b2 * b2_bounding_square) / b1_bounding_square
        max_b1 = int(min(limit_from_budget, limit_from_packing))

        if max_b1 < 0:
            continue

        for b1 in range(max_b1 + 1):
            # Calculate objective and cost for the current (b1, b2) pair
            current_objective = b1 * (b1_radius**2) + b2 * (b2_radius**2)
            current_cost = b1 * b1_cost + b2 * b2_cost
            
            # Check if this solution is better than the best one found so far
            if current_objective > max_objective_value:
                max_objective_value = current_objective
                best_b1 = b1
                best_b2 = b2
                highest_cost_for_max_objective = current_cost
            elif current_objective == max_objective_value:
                # Tie-breaker: prefer the solution that costs more (uses budget more fully)
                if current_cost > highest_cost_for_max_objective:
                    best_b1 = b1
                    best_b2 = b2
                    highest_cost_for_max_objective = current_cost

    # 3. Calculate and print the final results
    total_coverage = best_b1 * b1_area + best_b2 * b2_area
    coverage_ratio = total_coverage / city_area
    coverage_percentage = round(coverage_ratio * 100)
    final_cost = best_b1 * b1_cost + best_b2 * b2_cost

    print("--- Optimization Problem Formulation ---")
    print(f"Objective: Maximize Coverage = b1 * (\u03C0 * {b1_radius}\u00b2) + b2 * (\u03C0 * {b2_radius}\u00b2)")
    print("\n--- Constraints ---")
    print(f"Budget:    b1 * {b1_cost} + b2 * {b2_cost} \u2264 {total_budget}")
    print(f"Packing:   b1 * {int(b1_bounding_square)} + b2 * {int(b2_bounding_square)} \u2264 {city_area}")
    
    print("\n--- Optimal Solution ---")
    print(f"The optimal number of B1 towers (b1) is: {best_b1}")
    print(f"The optimal number of B2 towers (b2) is: {best_b2}")

    print("\n--- Resulting Equations ---")
    print(f"Total Cost:     {best_b1} * {b1_cost} + {best_b2} * {b2_cost} = {int(final_cost)} USD")
    print(f"Total Coverage: {best_b1} * \u03C0*{b1_radius}\u00b2 + {best_b2} * \u03C0*{b2_radius}\u00b2 = {total_coverage:.2f} km\u00b2")
    print(f"Coverage Ratio: ({total_coverage:.2f} km\u00b2 / {city_area} km\u00b2) * 100% = {coverage_percentage}%")
    
    # This is the final formatted answer for the system, not printed to the user.
    # The format is b1;b2;c
    # final_answer_string = f"{best_b1};{best_b2};{coverage_percentage}"


if __name__ == '__main__':
    solve_wifi_optimization()