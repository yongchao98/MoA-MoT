import math

def solve_tower_placement():
    """
    Solves the optimal tower placement problem to maximize coverage.
    """
    # Problem parameters
    city_width = 12  # km
    city_height = 11 # km
    budget = 45000   # usd

    cost_b1 = 1500
    radius_b1 = 1
    
    cost_b2 = 5000
    radius_b2 = 2

    # Calculate areas
    area_b1 = math.pi * (radius_b1 ** 2)
    area_b2 = math.pi * (radius_b2 ** 2)
    total_city_area = city_width * city_height

    # Variables to store the best solution found
    best_n1 = 0
    best_n2 = 0
    max_coverage = 0

    # Determine the maximum possible number of each tower type based on budget alone
    max_n1_possible = budget // cost_b1
    max_n2_possible = budget // cost_b2

    # Iterate through all possible combinations of n1 and n2
    for n2 in range(max_n2_possible + 1):
        for n1 in range(max_n1_possible + 1):
            # Calculate current cost and coverage for the combination (n1, n2)
            current_cost = n1 * cost_b1 + n2 * cost_b2
            current_coverage = n1 * area_b1 + n2 * area_b2

            # Check if the current combination is valid
            # 1. Is it within the budget?
            # 2. Is it physically possible to place (total area constraint)?
            if current_cost <= budget and current_coverage <= total_city_area:
                # If valid, check if this solution is better than the best one found so far
                if current_coverage > max_coverage:
                    max_coverage = current_coverage
                    best_n1 = n1
                    best_n2 = n2

    # Calculate final results for the best solution
    final_cost = best_n1 * cost_b1 + best_n2 * cost_b2
    coverage_ratio = max_coverage / total_city_area
    coverage_percentage = round(coverage_ratio * 100)

    # --- Outputting the formulation and solution details ---
    print("Optimization Problem Formulation:")
    print(f"Maximize Coverage = n1 * {area_b1:.2f} + n2 * {area_b2:.2f}")
    print("Subject to:")
    print(f"1. Budget: n1 * {cost_b1} + n2 * {cost_b2} <= {budget}")
    print(f"2. Area:   n1 * {area_b1:.2f} + n2 * {area_b2:.2f} <= {total_city_area}")
    print("\n--- Optimal Solution Found ---")
    print(f"Optimal number of B1 towers (n1): {best_n1}")
    print(f"Optimal number of B2 towers (n2): {best_n2}")

    print("\n--- Final Equation Values ---")
    print(f"Maximized Coverage = {best_n1} * {area_b1:.2f} + {best_n2} * {area_b2:.2f} = {max_coverage:.2f} km^2")
    print(f"Total Cost = {best_n1} * {cost_b1} + {best_n2} * {cost_b2} = {final_cost} USD")
    print(f"Coverage Ratio = ({max_coverage:.2f} / {total_city_area}) * 100 = {coverage_percentage}%")

    # Final answer in the required format
    print("\n--- Formatted Answer ---")
    final_answer = f"{best_n1};{best_n2};{coverage_percentage}"
    print(final_answer)


# Execute the function to solve the problem
solve_tower_placement()