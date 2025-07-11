import math

def solve_tower_optimization():
    """
    Solves the tower placement optimization problem based on a simplified model.
    """
    # --- Problem Parameters ---
    city_area = 12 * 11
    total_budget = 45000

    # B1 Tower specifications
    b1_cost = 1500
    b1_area = math.pi * (1**2)

    # B2 Tower specifications
    b2_cost = 5000
    b2_area = math.pi * (2**2)

    # --- Search for the Optimal Solution ---
    best_b1 = 0
    best_b2 = 0
    max_coverage = 0

    # Determine the maximum possible number of towers to limit the search space
    # Max b2 is limited by budget (45000 / 5000 = 9)
    max_b2_possible = int(total_budget / b2_cost)
    # Max b1 is limited by budget (45000 / 1500 = 30)
    max_b1_possible = int(total_budget / b1_cost)
    
    # Iterate through all feasible combinations of b1 and b2
    for b2 in range(max_b2_possible + 1):
        for b1 in range(max_b1_possible + 1):
            
            # Calculate total cost and total area for the current combination
            current_cost = b1 * b1_cost + b2 * b2_cost
            current_area = b1 * b1_area + b2 * b2_area

            # Check if the current combination satisfies both constraints
            if current_cost <= total_budget and current_area <= city_area:
                # If it's better than the best solution found so far, update it
                if current_area > max_coverage:
                    max_coverage = current_area
                    best_b1 = b1
                    best_b2 = b2

    # --- Calculate and Display Results ---
    final_cost = best_b1 * b1_cost + best_b2 * b2_cost
    coverage_ratio = (max_coverage / city_area) * 100
    rounded_coverage_ratio = round(coverage_ratio)
    
    print("Yes, the problem can be formulated for an optimization package, with a simplifying assumption for the placement constraint.")
    print("\n--- Simplified Problem Formulation ---")
    print("Let b1 = number of B1 towers, b2 = number of B2 towers.")
    print("\nObjective: Maximize Coverage")
    print(f"Maximize: b1 * {b1_area:.4f} + b2 * {b2_area:.4f}")

    print("\nConstraints:")
    print(f"1. Budget: b1 * {b1_cost} + b2 * {b2_cost} <= {total_budget}")
    print(f"2. Area:   b1 * {b1_area:.4f} + b2 * {b2_area:.4f} <= {city_area}")
    
    print("\n--- Solution ---")
    print(f"The optimal combination found is {best_b1} B1 towers and {best_b2} B2 towers.")

    print("\nVerification of the solution against the constraints:")
    # Per the instructions, showing each number in the final equation
    print(f"Final Cost: {best_b1} * {b1_cost} + {best_b2} * {b2_cost} = {final_cost} (<= {total_budget})")
    print(f"Final Area: {best_b1} * {b1_area:.4f} + {best_b2} * {b2_area:.4f} = {max_coverage:.4f} (<= {city_area})")

    print(f"\nThis results in a coverage ratio 'c' of {coverage_ratio:.2f}%, which is rounded to {rounded_coverage_ratio}%.")
    
    # Final output in the requested format b1;b2;c
    print("\n--- Final Answer ---")
    print(f"{best_b1};{best_b2};{rounded_coverage_ratio}")


# Execute the function to solve the problem
solve_tower_optimization()
<<<0;9;86>>>