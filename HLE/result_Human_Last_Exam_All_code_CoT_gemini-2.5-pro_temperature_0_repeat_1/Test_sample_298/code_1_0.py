import math

def solve_cohomology():
    """
    Calculates the cohomology groups of the moduli space M(7) and prints them.
    """
    k = 7
    # The moduli space M(k) is homotopy equivalent to the (k-1)-torus T^(k-1).
    n = k - 1

    cohomology_groups = []
    # The p-th cohomology group of T^n is Z^rank where rank is the binomial coefficient C(n, p).
    # The largest non-zero group is for p=n.
    for p in range(n + 1):
        # Calculate the rank of the cohomology group H^p
        rank = math.comb(n, p)
        
        # Format the group string according to the specified notation.
        if rank == 0:
            group_str = "0"
        elif rank == 1:
            group_str = "Z"
        else:
            group_str = f"Z^{rank}"
        cohomology_groups.append(group_str)
        
    # Format the final list as a string.
    result_str = "[" + ", ".join(cohomology_groups) + "]"
    print(result_str)

solve_cohomology()