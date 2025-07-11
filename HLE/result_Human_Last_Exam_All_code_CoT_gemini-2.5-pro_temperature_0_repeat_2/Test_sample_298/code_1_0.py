import math

def get_cohomology_group_str(rank):
    """Formats the cohomology group string based on its rank."""
    if rank == 0:
        return "0"
    elif rank == 1:
        return "Z"
    else:
        return f"Z^{rank}"

def solve_cohomology_m7():
    """
    Calculates and prints the list of cohomology groups for M(7).
    """
    k = 7
    cohomology_groups = []
    # The top non-zero cohomology group is H^k for a k-torus.
    for i in range(k + 1):
        # Calculate the rank using the binomial coefficient C(k, i)
        rank = math.comb(k, i)
        group_str = get_cohomology_group_str(rank)
        cohomology_groups.append(group_str)
    
    # Print the final list in the specified format
    print(f"[{', '.join(cohomology_groups)}]")

solve_cohomology_m7()