import math

def solve_tower_optimization():
    """
    Solves the tower placement optimization problem to maximize coverage.
    """
    # Define problem constants
    city_area = 12 * 11
    budget = 45000
    
    cost_b1 = 1500
    cost_b2 = 5000
    
    # Initialize variables to store the best solution found so far
    best_b1 = 0
    best_b2 = 0
    # The objective is to maximize b1 + 4*b2
    max_objective_value = -1

    # We can determine the maximum possible number of B2 towers from the budget
    max_b2 = math.floor(budget / cost_b2) # 45000 / 5000 = 9

    # Iterate through all possible numbers of B2 towers
    for b2 in range(max_b2 + 1):
        # For each number of b2, find the maximum number of b1 towers possible
        
        # 1. Calculate max b1 based on the remaining budget
        remaining_budget = budget - b2 * cost_b2
        max_b1_from_budget = math.floor(remaining_budget / cost_b1)
        
        # 2. Calculate max b1 based on the area constraint (b1 + 4*b2 <= 42)
        # Ensure the result is not negative before calculating
        if 42 - 4 * b2 >= 0:
            max_b1_from_area = 42 - 4 * b2
        else:
            continue # This combination of b2 is not possible

        # The number of b1 towers must satisfy both constraints
        current_b1 = min(max_b1_from_budget, max_b1_from_area)

        # Calculate the value of the objective function for this combination
        current_objective_value = current_b1 + 4 * b2
        
        # If this combination gives a better result, update the best solution
        if current_objective_value > max_objective_value:
            max_objective_value = current_objective_value
            best_b1 = current_b1
            best_b2 = b2
            
    # Once the loop is finished, we have the optimal number of towers
    b1_final = best_b1
    b2_final = best_b2

    # Calculate the final coverage area and ratio
    max_coverage_area = math.pi * max_objective_value
    coverage_ratio = (max_coverage_area / city_area) * 100
    
    # Rounded to the nearest percentage
    c_final = round(coverage_ratio)

    # Print the result in the format b1;b2;c
    # The final equation's numbers are the components of the answer itself.
    print(f"{b1_final};{b2_final};{c_final}")

solve_tower_optimization()