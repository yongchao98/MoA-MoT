import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # Define problem constants
    city_width = 12
    city_length = 11
    b1_radius = 1
    b2_radius = 2
    b1_cost = 1500
    b2_cost = 5000
    total_budget = 45000

    # Calculate areas
    city_area = city_width * city_length
    b1_area = math.pi * (b1_radius ** 2)
    b2_area = math.pi * (b2_radius ** 2)

    # Initialize variables to store the best solution found
    best_b1 = 0
    best_b2 = 0
    max_covered_area = 0.0

    # Determine the maximum possible number of each tower based on budget alone
    # This sets the search boundaries.
    max_b2_possible = total_budget // b2_cost

    # Iterate through all possible numbers of towers to find the optimal combination
    for b2 in range(max_b2_possible + 1):
        # For a given number of b2 towers, find the maximum b1 towers possible
        remaining_budget = total_budget - (b2 * b2_cost)
        max_b1_possible = remaining_budget // b1_cost
        
        for b1 in range(max_b1_possible + 1):
            # Calculate the total cost and total area for the current combination
            current_cost = (b1 * b1_cost) + (b2 * b2_cost)
            current_area = (b1 * b1_area) + (b2 * b2_area)

            # Check if the combination is valid (within budget and non-overlapping area)
            if current_cost <= total_budget and current_area <= city_area:
                # If valid, check if this combination provides better coverage
                if current_area > max_covered_area:
                    max_covered_area = current_area
                    best_b1 = b1
                    best_b2 = b2

    # Calculate the final coverage ratio as a percentage
    coverage_ratio = (max_covered_area / city_area) * 100
    # Round to the nearest whole percentage
    coverage_percentage = round(coverage_ratio)

    # Print the result in the specified format: b1;b2;c
    print(f"{best_b1};{best_b2};{coverage_percentage}")

# Execute the function to find and print the solution
solve_tower_placement()