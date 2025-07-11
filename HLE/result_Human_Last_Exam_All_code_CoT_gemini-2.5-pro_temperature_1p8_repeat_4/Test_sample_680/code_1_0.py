import math

def solve_tower_placement():
    """
    Solves the Wi-Fi tower placement optimization problem.
    """
    # Define problem parameters
    city_width = 12
    city_height = 11
    total_area = city_width * city_height

    b1_radius = 1
    b2_radius = 2

    # Initialize variables to store the best solution found
    best_n1 = 0
    best_n2 = 0
    max_objective_value = -1

    # Determine the maximum possible number of B2 towers based on constraints
    # From budget: 10*n2 <= 90 => n2 <= 9
    # From packing: 4*n2 <= 33 => n2 <= 8
    max_n2 = 8

    # Iterate through all possible numbers of B2 towers
    for n2 in range(max_n2 + 1):
        # For a given n2, find the maximum possible n1
        
        # From budget constraint: 3*n1 <= 90 - 10*n2
        n1_limit_budget = (90 - 10 * n2) / 3
        
        # From packing constraint: n1 <= 33 - 4*n2
        n1_limit_packing = 33 - 4 * n2

        # The number of n1 must be an integer and satisfy both constraints
        if n1_limit_budget < 0 or n1_limit_packing < 0:
            continue
            
        n1 = min(int(n1_limit_budget), int(n1_limit_packing))

        # Calculate the value of the objective function (proportional to coverage)
        current_objective_value = n1 + 4 * n2

        # If this solution is better or equal, update the best solution
        # (We use >= to find the solution with the most B2 towers in case of a tie)
        if current_objective_value >= max_objective_value:
            max_objective_value = current_objective_value
            best_n1 = n1
            best_n2 = n2

    # Calculate the final coverage ratio for the optimal solution
    coverage_area = math.pi * (best_n1 * b1_radius**2 + best_n2 * b2_radius**2)
    coverage_ratio_percent = round((coverage_area / total_area) * 100)

    # Output each number in the final equation
    print(f"Number of B1 towers (b1): {best_n1}")
    print(f"Number of B2 towers (b2): {best_n2}")
    print(f"Coverage Ratio (c): {coverage_ratio_percent}%")
    print("\nFinal formatted answer:")
    # Print the final answer in the requested format "b1;b2;c"
    print(f"{best_n1};{best_n2};{coverage_ratio_percent}")

solve_tower_placement()
<<<1;8;79>>>