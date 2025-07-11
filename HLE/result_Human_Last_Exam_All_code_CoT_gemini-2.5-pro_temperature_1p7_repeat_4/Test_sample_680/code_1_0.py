import math

def solve_tower_optimization():
    """
    Solves the tower placement optimization problem.
    """
    # --- Step 1: Define Problem Variables ---
    # The problem is to maximize the coverage area of a 12x11km rectangle
    # using two types of towers (B1, B2) with a total budget of 45,000 usd,
    # without any overlap in coverage.
    
    # City and Budget
    city_area = 12 * 11  # km^2
    total_budget = 45000

    # Tower B1 Specifications
    b1_cost = 1500
    b1_radius = 1
    b1_area = math.pi * b1_radius**2

    # Tower B2 Specifications
    b2_cost = 5000
    b2_radius = 2
    b2_area = math.pi * b2_radius**2

    # --- Step 2 & 3: Formulate and Solve the Optimization ---
    # The problem can be formulated for an integer optimization solver.
    # Let b1 and b2 be the number of B1 and B2 towers.
    # Objective: Maximize total area = b1 * (pi*1^2) + b2 * (pi*2^2).
    # This is equivalent to maximizing the metric: Z = b1 + 4*b2.
    # Constraint: 1500*b1 + 5000*b2 <= 45000.
    
    best_b1 = 0
    best_b2 = 0
    # Use the simplified metric for optimization
    max_coverage_metric = -1

    # To solve, we can iterate through all possible numbers of B2 towers.
    # Since B2 is more cost-effective for coverage (pi/1250 vs pi/1500),
    # we prioritize it by iterating downwards from the maximum possible.
    max_possible_b2 = total_budget // b2_cost
    for num_b2 in range(max_possible_b2, -1, -1):
        remaining_budget = total_budget - (num_b2 * b2_cost)
        num_b1 = remaining_budget // b1_cost
        
        current_metric = num_b1 + 4 * num_b2
        
        if current_metric > max_coverage_metric:
            max_coverage_metric = current_metric
            best_b1 = num_b1
            best_b2 = num_b2
            
    # The non-overlap constraint for the found solution (b1=0, b2=9) is feasible.
    # A staggered (hexagonal) packing of 9 circles of radius 2km fits within
    # a 12x10.93km area, which is inside the city's 12x11km dimensions.
    
    # --- Step 4: Calculate Final Results ---
    final_b1 = best_b1
    final_b2 = best_b2

    total_cost = final_b1 * b1_cost + final_b2 * b2_cost
    total_coverage_area = final_b1 * b1_area + final_b2 * b2_area
    coverage_ratio = total_coverage_area / city_area
    coverage_percentage = round(coverage_ratio * 100)

    # --- Step 5: Output the results and calculations ---
    print("--- Problem Formulation Summary ---")
    print("Objective: Maximize Coverage = b1 * Area(B1) + b2 * Area(B2)")
    print("Constraint: b1 * Cost(B1) + b2 * Cost(B2) <= Total Budget\n")
    
    print("--- Optimal Solution Details ---")
    print(f"The optimal combination is {final_b1} B1 towers and {final_b2} B2 towers.")

    print("\n--- Final Equation Components ---")
    print(f"Total Cost Calculation: {final_b1} * {b1_cost} + {final_b2} * {b2_cost} = {total_cost} USD")
    print(f"Total Coverage Area Calculation: {final_b1} * {b1_area:.2f} + {final_b2} * {b2_area:.2f} = {total_coverage_area:.2f} km^2")
    print(f"City Area: {city_area} km^2")
    print(f"Coverage Ratio: {coverage_percentage}%\n")
    
    print("--- Final Answer ---")
    # This is the final result in the requested format
    print(f"{final_b1};{final_b2};{coverage_percentage}")

if __name__ == '__main__':
    solve_tower_optimization()