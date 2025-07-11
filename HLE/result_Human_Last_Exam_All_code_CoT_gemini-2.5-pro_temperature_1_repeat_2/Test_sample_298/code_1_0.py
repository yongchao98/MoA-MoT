import math

def solve_cohomology():
    """
    This function calculates and prints the list of cohomology groups for M(7).
    """
    k = 7
    
    # The list of cohomology groups as strings
    cohomology_groups = []
    
    # H^0 is always Z for a connected space
    cohomology_groups.append("Z")
    
    # H^1 is Z based on the formula
    cohomology_groups.append("Z")
    
    # Calculate H^i for i >= 2
    n = k - 1
    # The highest non-trivial group is H^(k-1)
    for i in range(2, k):
        # The rank of H^i is the binomial coefficient C(k-1, i)
        rank = math.comb(n, i)
        
        if rank == 0:
            cohomology_groups.append("0")
        elif rank == 1:
            cohomology_groups.append("Z")
        else:
            cohomology_groups.append(f"Z^{rank}")
            
    # Format the final output string
    result_string = f"[{', '.join(cohomology_groups)}]"
    print(result_string)

solve_cohomology()