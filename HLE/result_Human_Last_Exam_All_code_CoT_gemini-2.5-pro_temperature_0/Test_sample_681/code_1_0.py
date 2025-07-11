import math

def solve_tower_optimization():
    """
    Calculates the optimal number of B1 and B2 towers to minimize cost
    while meeting the coverage requirements.
    """
    # Tower parameters
    cost_b1 = 1500
    cost_b2 = 4000
    coverage_unit_b1 = 1**2  # Represents t_i^2 for a B1 tower
    coverage_unit_b2 = 2**2  # Represents t_i^2 for a B2 tower
    required_coverage_units = 34

    min_cost = float('inf')
    best_b1 = -1
    best_b2 = -1

    # We can determine a reasonable upper bound for the number of B2 towers to check.
    # If we only used B1 towers, the cost would be 34 * 1500 = 51000.
    # The number of B2 towers will not exceed 51000 / 4000 = 12.75.
    # So, we can check for n2 from 0 to 13.
    max_n2_to_check = math.ceil((required_coverage_units * cost_b1) / cost_b2)

    for n2 in range(max_n2_to_check + 1):
        # For a given number of B2 towers (n2), find the minimum
        # number of B1 towers (n1) to meet the coverage constraint.
        remaining_coverage = required_coverage_units - n2 * coverage_unit_b2
        
        if remaining_coverage > 0:
            # n1 * coverage_unit_b1 >= remaining_coverage
            n1 = math.ceil(remaining_coverage / coverage_unit_b1)
        else:
            n1 = 0

        # Calculate the total cost for this combination of towers
        current_cost = n1 * cost_b1 + n2 * cost_b2

        # If this cost is the new minimum, store this solution
        if current_cost < min_cost:
            min_cost = current_cost
            best_b1 = n1
            best_b2 = n2

    # The problem formulation is correct, and this code finds the optimal solution.
    # The final output includes the number of each tower and the cost calculation.
    print(f"The optimal solution requires {best_b1} B1 towers and {best_b2} B2 towers.")
    print(f"The minimized cost is {best_b1} * {cost_b1} + {best_b2} * {cost_b2} = {min_cost}")
    print(f"<<<{best_b1};{best_b2};{min_cost}>>>")

solve_tower_optimization()