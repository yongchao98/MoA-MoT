import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem by iterating through all
    feasible options and finding the one that maximizes coverage.
    """
    # City and Budget parameters
    city_width = 12  # km
    city_length = 11  # km
    city_area = city_width * city_length
    total_budget = 45000

    # Tower B1 parameters
    cost1 = 1500
    radius1 = 1
    coverage1 = math.pi * radius1**2
    # Bounding box area for non-overlapping approximation
    footprint1 = (2 * radius1)**2

    # Tower B2 parameters
    cost2 = 5000
    radius2 = 2
    coverage2 = math.pi * radius2**2
    # Bounding box area for non-overlapping approximation
    footprint2 = (2 * radius2)**2

    # Variables to store the optimal solution found
    best_n1 = 0
    best_n2 = 0
    max_coverage = 0

    # Determine the maximum possible number of each tower type based on budget
    max_b1_possible = total_budget // cost1
    max_b2_possible = total_budget // cost2

    # Iterate through all combinations of n1 (B1 towers) and n2 (B2 towers)
    for n2 in range(max_b2_possible + 1):
        for n1 in range(max_b1_possible + 1):
            # 1. Check budget constraint
            current_cost = n1 * cost1 + n2 * cost2
            if current_cost > total_budget:
                # If budget is exceeded, further increase in n1 is not possible for this n2
                break

            # 2. Check placement constraint (approximated)
            current_footprint = n1 * footprint1 + n2 * footprint2
            if current_footprint > city_area:
                continue

            # Calculate total coverage for this valid combination
            current_coverage = n1 * coverage1 + n2 * coverage2

            # If this combination yields better coverage, update the best solution
            if current_coverage > max_coverage:
                max_coverage = current_coverage
                best_n1 = n1
                best_n2 = n2

    # Calculate the coverage ratio as a percentage for the best solution
    coverage_ratio = (max_coverage / city_area) * 100
    # Round the ratio to the nearest whole percentage point
    rounded_coverage_percentage = round(coverage_ratio)

    # Print the final result in the format b1;b2;c
    # This output provides the number of B1 towers, B2 towers, and the coverage percentage.
    print(f"{best_n1};{best_n2};{rounded_coverage_percentage}")

if __name__ == "__main__":
    solve_tower_placement()