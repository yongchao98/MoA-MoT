import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # Define problem parameters
    city_area = 12 * 11  # km^2
    budget = 45000       # usd

    # B1 tower properties
    cost_b1 = 1500
    radius_b1 = 1
    area_b1 = math.pi * (radius_b1 ** 2)

    # B2 tower properties
    cost_b2 = 5000
    radius_b2 = 2
    area_b2 = math.pi * (radius_b2 ** 2)

    best_n1 = 0
    best_n2 = 0
    max_coverage_factor = -1

    # Determine the search range based on budget
    max_n1 = budget // cost_b1
    max_n2 = budget // cost_b2

    # Iterate through all possible combinations of n1 and n2
    for n1 in range(max_n1 + 1):
        for n2 in range(max_n2 + 1):
            
            # Check budget constraint
            total_cost = (n1 * cost_b1) + (n2 * cost_b2)
            if total_cost > budget:
                continue

            # Check non-overlapping area constraint
            total_covered_area = (n1 * area_b1) + (n2 * area_b2)
            if total_covered_area > city_area:
                continue
            
            # The objective is to maximize π * (n1 + 4*n2), which is equivalent to maximizing (n1 + 4*n2)
            current_coverage_factor = n1 + 4 * n2
            
            if current_coverage_factor > max_coverage_factor:
                max_coverage_factor = current_coverage_factor
                best_n1 = n1
                best_n2 = n2

    # The optimal number of towers
    b1 = best_n1
    b2 = best_n2
    
    # Calculate the final coverage percentage
    optimal_coverage_area = (b1 * area_b1) + (b2 * area_b2)
    coverage_ratio = optimal_coverage_area / city_area
    c = round(coverage_ratio * 100)

    # Print each number in the final equation as requested
    print("--- Optimal Solution ---")
    print(f"Number of B1 towers (b1): {b1}")
    print(f"Number of B2 towers (b2): {b2}")
    print(f"Resulting coverage percentage (c): {c}")
    
    # Final answer in the required format
    final_answer = f"{b1};{b2};{c}"
    print(f"\nFinal answer string in b1;b2;c format:")
    print(f"<<<{final_answer}>>>")

# Run the solver
solve_tower_placement()