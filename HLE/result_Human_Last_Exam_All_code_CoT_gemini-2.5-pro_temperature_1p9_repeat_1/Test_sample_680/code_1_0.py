import math

def solve_tower_placement():
    """
    Solves the tower placement optimization problem to maximize coverage
    within a given budget and area constraints.
    """
    # Define problem constants
    city_area = 12 * 11
    
    # Tower B1 properties
    b1_radius = 1
    b1_cost = 1500
    b1_area = math.pi * b1_radius**2

    # Tower B2 properties
    b2_radius = 2
    b2_cost = 5000
    b2_area = math.pi * b2_radius**2

    total_budget = 45000

    # Variables to store the optimal solution found
    best_b1 = 0
    best_b2 = 0
    max_coverage = 0.0

    # Determine the maximum possible number of B2 towers based on budget
    max_b2_towers = total_budget // b2_cost

    # Iterate through all possible numbers of B2 towers (from 0 to max_b2_towers)
    for b2 in range(max_b2_towers + 1):
        # Calculate the remaining budget after placing b2 towers
        remaining_budget = total_budget - (b2 * b2_cost)
        
        # Calculate the maximum number of B1 towers that can be afforded
        b1 = remaining_budget // b1_cost
        
        # Calculate the total coverage for the current combination of towers
        current_coverage = (b1 * b1_area) + (b2 * b2_area)
        
        # The simplified non-overlapping constraint: total circle area must be <= city area
        if current_coverage <= city_area:
            # If this combination gives a better coverage, we update our best solution
            if current_coverage > max_coverage:
                max_coverage = current_coverage
                best_b1 = b1
                best_b2 = b2
    
    # Calculate the final coverage ratio as a percentage
    coverage_ratio = (max_coverage / city_area) * 100
    rounded_coverage_ratio = round(coverage_ratio)
    
    # --- Output the results ---
    print(f"Based on the optimization, the recommended number of towers are:")
    print(f"- B1 Towers (b1): {best_b1}")
    print(f"- B2 Towers (b2): {best_b2}")
    
    print("\nThis combination results in:")
    total_cost = best_b1 * b1_cost + best_b2 * b2_cost
    print(f"- Total Cost: ${total_cost:,.0f} (Budget: ${total_budget:,.0f})")
    print(f"- Total Coverage Area: {max_coverage:.2f} km² (City Area: {city_area} km²)")
    
    print("\nThe final answer in the format b1;b2;c is:")
    # Final equation format as requested
    print(f"{best_b1};{best_b2};{int(rounded_coverage_ratio)}")

# Run the solver
solve_tower_placement()