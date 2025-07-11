import math

def solve_tower_placement():
    """
    Formulates and solves the tower placement optimization problem.
    """
    # Tower and city specifications
    radius_b1 = 1
    cost_b1 = 1500
    area_factor_b1 = radius_b1**2

    radius_b2 = 2
    cost_b2 = 5000
    area_factor_b2 = radius_b2**2

    city_area = 12 * 11
    total_budget = 45000

    # --- Problem Formulation Output ---
    print("Step 1: Formulating the Optimization Problem")
    print("Let b1 be the number of B1 towers and b2 be the number of B2 towers.")
    
    print("\nObjective Function (Maximize Total Area):")
    print(f"Maximize: b1 * pi * {radius_b1}^2 + b2 * pi * {radius_b2}^2")
    print(f"This is equivalent to maximizing Z = {area_factor_b1}*b1 + {area_factor_b2}*b2")

    print("\nConstraints:")
    # Budget constraint
    simplified_cost_b1 = int(cost_b1 / 500)
    simplified_cost_b2 = int(cost_b2 / 500)
    simplified_budget = int(total_budget / 500)
    print(f"1. Budget: {cost_b1}*b1 + {cost_b2}*b2 <= {total_budget}")
    print(f"   Simplified: {simplified_cost_b1}*b1 + {simplified_cost_b2}*b2 <= {simplified_budget}")

    # Area constraint
    area_constraint_rhs = int(city_area / math.pi)
    print(f"2. Non-overlapping Area: b1 * pi * {radius_b1}^2 + b2 * pi * {radius_b2}^2 <= {city_area}")
    print(f"   Simplified: {area_factor_b1}*b1 + {area_factor_b2}*b2 <= {area_constraint_rhs}")
    
    print("-" * 30)

    # --- Problem Solving ---
    print("\nStep 2: Solving the problem by iterating through feasible solutions")
    best_b1 = 0
    best_b2 = 0
    max_z = 0

    # Iterate through possible numbers of B2 towers
    # Max b2 is limited by budget (45000/5000=9) and area (42/4=10.5)
    max_b2_budget = total_budget // cost_b2
    for b2 in range(max_b2_budget + 1):
        # For each b2, find the max b1 possible under the constraints
        
        # From budget constraint
        remaining_budget = total_budget - b2 * cost_b2
        max_b1_budget = remaining_budget // cost_b1

        # From area constraint
        # b1 <= 42 - 4*b2
        remaining_area_factor = area_constraint_rhs - b2 * area_factor_b2
        if remaining_area_factor < 0:
            max_b1_area = -1 # Impossible
        else:
            max_b1_area = remaining_area_factor // area_factor_b1
        
        # The number of b1 towers must satisfy both constraints
        b1 = min(max_b1_budget, max_b1_area)

        if b1 < 0:
            continue

        # Calculate the objective function Z
        current_z = b1 * area_factor_b1 + b2 * area_factor_b2
        
        if current_z > max_z:
            max_z = current_z
            best_b1 = b1
            best_b2 = b2

    # --- Final Result ---
    # Calculate coverage ratio
    total_coverage_area = max_z * math.pi
    coverage_ratio = total_coverage_area / city_area
    coverage_percentage = round(coverage_ratio * 100)
    
    print("\nStep 3: Final Answer")
    print(f"The optimal solution is to use {best_b1} B1 towers and {best_b2} B2 towers.")
    print(f"This results in a coverage ratio of {coverage_percentage}%.")

    # Final output in the required format
    print("\nResult in b1;b2;c format:")
    print(f"{best_b1};{best_b2};{coverage_percentage}")

# Execute the function
solve_tower_placement()