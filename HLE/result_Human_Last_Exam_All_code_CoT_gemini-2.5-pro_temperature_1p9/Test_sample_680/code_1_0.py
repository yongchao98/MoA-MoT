import math

def solve_tower_optimization():
    """
    Formulates and solves the tower placement problem to maximize coverage
    within a given budget, based on a simplified Integer Linear Program.
    """
    
    # --- 1. Define Variables and Constants ---
    total_budget = 45000
    city_area_total = 12 * 11  # in km^2

    # Tower B1 parameters
    b1_cost = 1500
    b1_radius = 1
    b1_area = math.pi * b1_radius**2

    # Tower B2 parameters
    b2_cost = 5000
    b2_radius = 2
    b2_area = math.pi * b2_radius**2

    # --- 2 & 3. Initialize for ILP solution ---
    # We are looking for non-negative integers b1, b2 to:
    # Maximize: Z = b1 + 4*b2
    # Subject to: 1500*b1 + 5000*b2 <= 45000
    
    best_b1 = 0
    best_b2 = 0
    max_objective_value = -1

    # --- 4. Solve the ILP by iterating through possibilities ---
    # We iterate through the more expensive tower (B2) to reduce loops.
    max_b2_count = total_budget // b2_cost
    
    for num_b2 in range(max_b2_count + 1):
        # Calculate the budget remaining after purchasing B2 towers
        remaining_budget = total_budget - (num_b2 * b2_cost)
        
        # Calculate the maximum number of B1 towers we can afford
        num_b1 = remaining_budget // b1_cost
        
        # Calculate the value of our objective function
        current_objective_value = num_b1 + 4 * num_b2
        
        # If this combination gives a better result, we store it
        if current_objective_value > max_objective_value:
            max_objective_value = current_objective_value
            best_b1 = num_b1
            best_b2 = num_b2
            
    # --- 5. Calculate and Print the Final Answer ---
    
    # Calculate the total coverage area with the optimal number of towers
    final_coverage_area = (best_b1 * b1_area) + (best_b2 * b2_area)
    
    # Calculate the coverage ratio as a percentage
    coverage_ratio = (final_coverage_area / city_area_total) * 100
    
    # Round the coverage ratio to the nearest whole number for the output
    rounded_coverage_percentage = int(round(coverage_ratio))
    
    # Print the final result in the format b1;b2;c
    print(f"{best_b1};{best_b2};{rounded_coverage_percentage}")

# Execute the function
solve_tower_optimization()
<<<0;9;86>>>