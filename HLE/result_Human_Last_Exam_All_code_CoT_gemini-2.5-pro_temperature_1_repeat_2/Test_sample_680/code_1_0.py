import math

def solve_tower_optimization():
    """
    Solves the tower placement optimization problem to maximize coverage area
    within a given budget and area.
    """
    # --- Problem Definition ---
    # City dimensions
    city_width = 12  # km
    city_length = 11 # km
    city_area = city_width * city_length

    # Tower specifications
    b1_radius = 1.0
    b1_cost = 1500
    b1_area = math.pi * b1_radius**2

    b2_radius = 2.0
    b2_cost = 5000
    b2_area = math.pi * b2_radius**2

    # Budget constraint
    total_budget = 45000

    # --- Optimization ---
    # We want to find integers b1, b2 >= 0 that:
    # Maximize: b1 * b1_area + b2 * b2_area  (which is proportional to b1 + 4*b2)
    # Subject to: b1 * b1_cost + b2 * b2_cost <= total_budget

    best_b1 = 0
    best_b2 = 0
    max_objective_value = -1

    # Iterate through all possible numbers of the more expensive tower (B2)
    max_b2_possible = total_budget // b2_cost
    for b2 in range(max_b2_possible + 1):
        remaining_budget = total_budget - (b2 * b2_cost)
        
        # Calculate max number of B1 towers with the remaining budget
        b1 = remaining_budget // b1_cost
        
        # This is the value we want to maximize
        objective_value = b1 * 1 + b2 * 4
        
        if objective_value > max_objective_value:
            max_objective_value = objective_value
            best_b1 = b1
            best_b2 = b2

    # --- Calculate Final Results ---
    b1 = best_b1
    b2 = best_b2
    
    total_cost = b1 * b1_cost + b2 * b2_cost
    total_coverage = b1 * b1_area + b2 * b2_area
    coverage_ratio = total_coverage / city_area
    coverage_percentage = round(coverage_ratio * 100)
    
    # --- Output the results ---
    print("--- Optimization Formulation & Solution ---")
    print(f"Objective: Maximize Coverage = b1 * {b1_area:.2f} + b2 * {b2_area:.2f}")
    print(f"Equivalent Objective: Maximize Z = b1 * 1 + b2 * 4")
    print(f"Constraint: b1 * {b1_cost} + b2 * {b2_cost} <= {total_budget}\n")

    print(f"Optimal number of B1 towers (b1): {b1}")
    print(f"Optimal number of B2 towers (b2): {b2}\n")

    print("--- Final Equation Numbers ---")
    print(f"Maximized Objective Equation: Z = {b1}*1 + {b2}*4 = {b1*1 + b2*4}")
    print(f"Budget Equation: {b1}*{b1_cost} + {b2}*{b2_cost} = {total_cost} <= {total_budget}")
    print(f"Coverage Equation: {b1} * {b1_area:.2f} + {b2} * {b2_area:.2f} = {total_coverage:.2f} km^2")
    print(f"City Area: {city_area} km^2")
    print(f"Coverage Ratio: ({total_coverage:.2f} / {city_area:.2f}) * 100 = {coverage_ratio*100:.2f}%")
    print(f"Final Rounded Coverage Percentage (c): {coverage_percentage}%\n")

    # The final answer in the required format
    final_answer = f"{b1};{b2};{coverage_percentage}"
    print(f"Final Answer String: {final_answer}")
    
    # Do not remove the following line, it is for the final answer.
    print(f"<<<{final_answer}>>>")

solve_tower_optimization()