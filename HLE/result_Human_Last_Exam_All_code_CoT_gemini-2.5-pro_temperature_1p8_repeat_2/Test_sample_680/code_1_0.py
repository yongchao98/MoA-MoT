import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem.
    """
    # --- 1. Define problem constants ---
    city_width = 12
    city_length = 11
    total_city_area = city_width * city_length
    total_budget = 45000

    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * (radius_b1**2)

    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * (radius_b2**2)

    # --- Print problem formulation ---
    print("This problem can be formulated as an integer optimization problem.")
    print("\nVariables:")
    print("  b1: number of B1 towers (radius 1km, cost 1,500)")
    print("  b2: number of B2 towers (radius 2km, cost 5,000)")

    print("\nObjective: Maximize Total Coverage Area")
    # Using unicode for pi
    print(f"  Maximize: b1 * \u03c0*(1^2) + b2 * \u03c0*(2^2)  = \u03c0 * (b1 + 4*b2)")

    print("\nConstraints:")
    print("1. Budget Constraint:")
    print(f"   {cost_b1}*b1 + {cost_b2}*b2 <= {total_budget}")
    print(f"   Simplified: {cost_b1 // 500}*b1 + {cost_b2 // 500}*b2 <= {total_budget // 500}")

    print("\n2. Non-overlapping Area Constraint (simplified):")
    print("   The sum of tower areas must be less than the total city area.")
    print(f"   \u03c0 * (b1 + 4*b2) <= {city_width} * {city_length}")
    print(f"   b1 + 4*b2 <= {total_city_area} / \u03c0")
    print(f"   b1 + 4*b2 <= {math.floor(total_city_area / math.pi)}")
    print("\n----------------------------------------------------")
    
    # --- 2. Initialize variables for the search ---
    max_objective_value = -1
    best_b1 = 0
    best_b2 = 0

    # --- 3. Brute-force search for the optimal solution ---
    # Determine the maximum possible number of towers based on budget alone
    max_b1_possible = total_budget // cost_b1
    max_b2_possible = total_budget // cost_b2

    for b2 in range(max_b2_possible + 1):
        for b1 in range(max_b1_possible + 1):
            # Check budget constraint (simplified version)
            if (3 * b1 + 10 * b2) > 90:
                continue

            # Check area constraint (simplified version)
            if (b1 + 4 * b2) > math.floor(total_city_area / math.pi):
                continue

            # This is a valid combination, calculate its objective value
            current_objective_value = b1 + 4 * b2
            
            # If this solution is better than the best one found so far, update
            if current_objective_value > max_objective_value:
                max_objective_value = current_objective_value
                best_b1 = b1
                best_b2 = b2

    # --- 4. Calculate final results ---
    max_coverage_area = best_b1 * area_b1 + best_b2 * area_b2
    coverage_ratio = (max_coverage_area / total_city_area) * 100
    coverage_percentage = round(coverage_ratio)
    
    # --- 5. Print the final answer ---
    print("\nOptimal Solution:")
    print(f"Number of B1 towers (b1): {best_b1}")
    print(f"Number of B2 towers (b2): {best_b2}")
    
    final_cost = best_b1 * cost_b1 + best_b2 * cost_b2
    print(f"Total cost: ${final_cost:,} (Budget: ${total_budget:,})")
    
    print(f"Maximized coverage area: {max_coverage_area:.2f} km^2")
    print(f"Coverage Ratio: {coverage_ratio:.2f}% of the {total_city_area} km^2 city area.")
    
    print("\nFinal Answer (format: b1;b2;coverage_percentage):")
    result_string = f"{best_b1};{best_b2};{coverage_percentage}"
    print(result_string)

# Run the solver
solve_tower_placement()