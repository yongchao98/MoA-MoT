import math

def solve_tower_optimization():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # Problem Parameters
    total_area = 12 * 11  # km^2
    budget = 45000  # usd

    # Tower B1 parameters
    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * radius_b1**2

    # Tower B2 parameters
    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * radius_b2**2
    
    # Simplified budget constraint: 3*b1 + 10*b2 <= 90
    
    best_b1 = 0
    best_b2 = 0
    max_coverage_metric = -1

    # Iterate through all possible numbers of B2 towers
    # The max number of b2 is floor(90 / 10) = 9
    max_b2_possible = int(90 / 10)
    for b2 in range(max_b2_possible, -1, -1):
        # For a given b2, calculate the max number of b1 towers allowed by the budget
        remaining_budget_units = 90 - 10 * b2
        b1 = int(remaining_budget_units / 3)
        
        # Calculate the coverage metric C = b1 + 4*b2
        coverage_metric = b1 + 4 * b2

        # Check the area constraint: b1 + 4*b2 <= 42
        # Note: coverage_metric is the left side of this inequality
        if coverage_metric * math.pi > total_area:
            # This combination is not physically possible, skip
            continue

        # If this combination gives better coverage, store it
        if coverage_metric > max_coverage_metric:
            max_coverage_metric = coverage_metric
            best_b1 = b1
            best_b2 = b2

    # Calculate the final coverage ratio with the best combination
    total_covered_area = (best_b1 * area_b1) + (best_b2 * area_b2)
    coverage_ratio = (total_covered_area / total_area) * 100
    
    # Round to the nearest percentage
    final_c = round(coverage_ratio)

    # Print the results as requested
    print(f"--- Optimization Results ---")
    print(f"Optimal number of B1 towers (b1): {best_b1}")
    print(f"Optimal number of B2 towers (b2): {best_b2}")
    print("\n--- Coverage Calculation ---")
    print(f"Total Covered Area = ({best_b1} * {area_b1:.2f}) + ({best_b2} * {area_b2:.2f}) = {total_covered_area:.2f} km^2")
    print(f"Total City Area = {total_area} km^2")
    print(f"Coverage Ratio (c) = ({total_covered_area:.2f} / {total_area}) * 100 = {coverage_ratio:.2f}%")
    print(f"Rounded Coverage Percentage: {final_c}%")
    
    # Final answer in the required format
    print("\n--- Final Answer Format ---")
    print(f"{best_b1};{best_b2};{final_c}")

solve_tower_optimization()