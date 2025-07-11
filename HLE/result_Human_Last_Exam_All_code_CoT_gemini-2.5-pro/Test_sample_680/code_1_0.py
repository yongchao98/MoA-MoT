import math

def solve_tower_optimization():
    """
    This function solves the tower placement optimization problem by formulating it as an
    integer linear program and finding the optimal number of towers.
    """
    # Define problem parameters
    b1_radius = 1
    b2_radius = 2
    b1_cost = 1500
    b2_cost = 5000
    total_budget = 45000
    city_width = 12
    city_length = 11

    # --- Problem Formulation ---
    # Objective: Maximize Coverage Area = b1*pi*r1^2 + b2*pi*r2^2
    # This is equivalent to maximizing the expression: Z = b1 + 4*b2

    # Constraints:
    # 1. Budget: b1_cost*b1 + b2_cost*b2 <= total_budget
    #    => 1500*b1 + 5000*b2 <= 45000
    #    => 3*b1 + 10*b2 <= 90
    #
    # 2. Placement (non-overlap using bounding squares):
    #    b1*(2*r1)^2 + b2*(2*r2)^2 <= city_width * city_length
    #    => 4*b1 + 16*b2 <= 132
    #    => b1 + 4*b2 <= 33

    # --- Solver ---
    # Iterate through all possible integer solutions to find the optimum.
    # From the placement constraint, 4*b2 <= 33, so b2_max is 8.
    b2_max = (city_width * city_length) // ((2 * b2_radius)**2)
    
    max_objective_value = -1
    optimal_solutions = []

    for b2 in range(b2_max + 1):
        # From budget constraint: 3*b1 <= 90 - 10*b2
        n1_from_budget = (90 - 10 * b2) // 3
        
        # From placement constraint: b1 <= 33 - 4*b2
        n1_from_placement = 33 - 4 * b2

        # The number of b1 towers is limited by the stricter of the two constraints
        # and must be non-negative.
        if n1_from_budget < 0 or n1_from_placement < 0:
            continue
        b1 = min(n1_from_budget, n1_from_placement)

        # Calculate the value of the objective function
        current_objective = b1 + 4 * b2

        if current_objective > max_objective_value:
            max_objective_value = current_objective
            optimal_solutions = [(b1, b2)]
        elif current_objective == max_objective_value:
            optimal_solutions.append((b1, b2))

    # --- Tie-Breaking ---
    # If multiple solutions give the same max coverage, choose the one with the highest cost.
    best_b1, best_b2 = -1, -1
    max_cost = -1
    for b1_sol, b2_sol in optimal_solutions:
        cost = b1_cost * b1_sol + b2_cost * b2_sol
        if cost > max_cost:
            max_cost = cost
            best_b1 = b1_sol
            best_b2 = b2_sol

    # --- Final Calculation ---
    total_coverage_area = best_b1 * math.pi * (b1_radius**2) + best_b2 * math.pi * (b2_radius**2)
    total_city_area = city_width * city_length
    coverage_ratio = total_coverage_area / total_city_area
    coverage_percentage = round(coverage_ratio * 100)

    # --- Output ---
    # Print the final result in the format b1;b2;c
    print(f"{best_b1};{best_b2};{coverage_percentage}")

solve_tower_optimization()
<<<13;5;79>>>