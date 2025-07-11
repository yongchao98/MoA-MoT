import math

def solve_tower_optimization():
    """
    Formulates and solves the tower placement optimization problem.
    """
    # Problem constants
    city_width = 12.0  # km
    city_height = 11.0 # km
    total_budget = 45000 # usd

    # Tower B1 properties
    radius1 = 1.0  # km
    cost1 = 1500   # usd
    area1 = math.pi * radius1**2

    # Tower B2 properties
    radius2 = 2.0  # km
    cost2 = 5000   # usd
    area2 = math.pi * radius2**2

    city_area = city_width * city_height

    # --- Problem Formulation Output ---
    print("Optimization Problem Formulation:")
    print("---------------------------------")
    print("Variables:")
    print("  b1: Number of B1 towers (radius 1km, cost 1500)")
    print("  b2: Number of B2 towers (radius 2km, cost 5000)")
    print("\nObjective: Maximize Total Coverage")
    print(f"  Maximize: b1 * ({area1:.2f}) + b2 * ({area2:.2f})")
    
    print("\nConstraints:")
    print(f"  1. Budget: b1 * {cost1} + b2 * {cost2} <= {total_budget}")
    print(f"  2. Area:   b1 * ({area1:.2f}) + b2 * ({area2:.2f}) <= {city_area}")
    print("---------------------------------\n")


    # --- Solving the problem by iterating through possibilities ---
    best_b1 = 0
    best_b2 = 0
    max_coverage_area = 0.0

    # Determine the maximum number of B2 towers we can afford
    max_b2_possible = int(total_budget / cost2)

    # Iterate through all possible numbers of B2 towers
    for b2 in range(max_b2_possible + 1):
        # For a given number of b2 towers, calculate remaining budget
        remaining_budget = total_budget - (b2 * cost2)
        
        # Calculate the max number of b1 towers we can afford with the rest
        b1 = int(remaining_budget / cost1)
        
        # Calculate current combination's total area and cost
        current_area = b1 * area1 + b2 * area2
        current_cost = b1 * cost1 + b2 * cost2
        
        # Check if this combination respects the area constraint
        if current_area <= city_area:
            # If it's the best coverage so far, store it
            if current_area > max_coverage_area:
                max_coverage_area = current_area
                best_b1 = b1
                best_b2 = b2

    # --- Final Result ---
    # Calculate the final coverage ratio
    coverage_ratio = max_coverage_area / city_area
    coverage_percentage = round(coverage_ratio * 100)

    # Print the details of the solution
    print("Optimal Solution Found:")
    print(f"  Number of B1 towers (b1): {best_b1}")
    print(f"  Number of B2 towers (b2): {best_b2}\n")
    
    print("Final Equation Values:")
    final_cost = best_b1 * cost1 + best_b2 * cost2
    print(f"  Total Cost = {best_b1} * {cost1} + {best_b2} * {cost2} = {final_cost} (<= {total_budget})")
    print(f"  Total Area = {best_b1} * {area1:.2f} + {best_b2} * {area2:.2f} = {max_coverage_area:.2f} (<= {city_area})")
    print(f"  Coverage Ratio = {max_coverage_area:.2f} / {city_area:.2f} = {coverage_ratio:.4f}")
    print(f"  Coverage Percentage (c) = {coverage_percentage}%\n")


    # Output the final answer in the required format
    final_answer = f"{best_b1};{best_b2};{coverage_percentage}"
    print("Formatted Answer:")
    print(final_answer)
    return final_answer

final_answer_string = solve_tower_optimization()
print(f"<<<{final_answer_string}>>>")