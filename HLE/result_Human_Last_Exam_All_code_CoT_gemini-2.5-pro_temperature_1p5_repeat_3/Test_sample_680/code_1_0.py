import math

def solve_tower_placement():
    """
    Solves the WiFi tower placement optimization problem.
    """
    # --- Problem Definition ---
    # City dimensions and total area
    city_width = 12  # km
    city_length = 11 # km
    city_area = city_width * city_length

    # Tower specifications
    b1_cost = 1500
    b1_radius = 1
    b1_area = math.pi * b1_radius**2

    b2_cost = 5000
    b2_radius = 2
    b2_area = math.pi * b2_radius**2

    # Total budget
    total_budget = 45000

    # --- Solving the Optimization Problem ---
    best_b1 = 0
    best_b2 = 0
    max_coverage = 0

    # Determine the search range for b2 towers.
    # The number of b2 towers cannot exceed the budget or the total area.
    max_b2_by_budget = int(total_budget / b2_cost)
    max_b2_by_area = int(city_area / b2_area)
    max_b2 = min(max_b2_by_budget, max_b2_by_area)

    # Iterate through all possible numbers of b2 towers
    for num_b2 in range(max_b2 + 1):
        # Calculate remaining budget and area after placing b2 towers
        remaining_budget = total_budget - num_b2 * b2_cost
        remaining_area = city_area - num_b2 * b2_area

        # If remaining budget or area is negative, this combination is invalid
        if remaining_budget < 0 or remaining_area < 0:
            continue

        # Calculate max possible b1 towers with remaining resources
        max_b1_by_budget = int(remaining_budget / b1_cost)
        max_b1_by_area = int(remaining_area / b1_area)

        # The actual number of b1 towers is limited by both budget and area
        num_b1 = min(max_b1_by_budget, max_b1_by_area)

        # Calculate the total coverage for this combination
        current_coverage = num_b1 * b1_area + num_b2 * b2_area

        # If this combination provides better coverage, store it as the best
        if current_coverage > max_coverage:
            max_coverage = current_coverage
            best_b1 = num_b1
            best_b2 = num_b2

    # --- Final Result Calculation ---
    # Calculate coverage ratio as a percentage
    coverage_ratio = (max_coverage / city_area) * 100
    # Round to the nearest whole percentage point
    coverage_percentage = int(round(coverage_ratio))

    # Print the result in the specified format
    print(f"{best_b1};{best_b2};{coverage_percentage}")

# Execute the function to find and print the solution
solve_tower_placement()