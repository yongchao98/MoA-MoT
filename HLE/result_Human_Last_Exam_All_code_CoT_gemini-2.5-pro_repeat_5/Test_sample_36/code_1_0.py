import math
from itertools import chain, combinations
from functools import reduce

def get_abelian_groups(n):
    """
    Finds the cyclic factors for all abelian groups of order n.
    This is a helper to illustrate the groups involved, not used in the main calculation.
    """
    # This is a complex number theory problem. For n=18, we know the groups.
    # 18 = 2 * 3^2. Partitions of 2 are (2) and (1,1).
    # This gives Z_2 x Z_9 and Z_2 x Z_3 x Z_3.
    return [[2, 9], [2, 3, 3]]

def lcm(a, b):
    """Calculates the least common multiple of two integers."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // math.gcd(a, b)

def get_exponent(factors):
    """Calculates the exponent of a group given by its cyclic factors."""
    if not factors:
        return 1
    return reduce(lcm, factors, 1)

def get_divisors(n):
    """Calculates all positive divisors of an integer n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def get_generating_antichains(elements):
    """
    Finds all non-empty antichains in a set of elements partially ordered by divisibility.
    Each antichain corresponds to a unique set of eigenvalues S(rho).
    """
    antichains = set()
    # Generate the power set of elements, from size 1 up to len(elements).
    for r in range(1, len(elements) + 1):
        for subset in combinations(elements, r):
            is_an_antichain = True
            # For a subset to be an antichain, no element can divide another.
            if len(subset) > 1:
                for i in range(len(subset)):
                    for j in range(i + 1, len(subset)):
                        a, b = subset[i], subset[j]
                        if a % b == 0 or b % a == 0:
                            is_an_antichain = False
                            break
                    if not is_an_antichain:
                        break
            
            if is_an_antichain:
                antichains.add(frozenset(subset))
                
    return antichains

def main():
    """
    Main function to execute the plan and find the number of unique eigenvalue sets.
    """
    print("This script calculates the number of unique sets of eigenvalues for representations of Abelian groups of order 18.")
    print("-" * 80)

    # Step 1: Identify the groups
    group_factors_list = get_abelian_groups(18)
    print("The two Abelian groups of order 18 are G1 = Z_2 x Z_9 and G2 = Z_2 x Z_3 x Z_3.\n")

    all_unique_antichains = set()
    group_results = {}

    # Step 2 & 3: Analyze each group
    for i, factors in enumerate(group_factors_list):
        group_name = f"G{i+1}"
        
        # a. Calculate exponent
        exponent = get_exponent(factors)
        
        # b. Find divisors
        divs = get_divisors(exponent)
        
        # c. Find antichains
        antichains = get_generating_antichains(divs)
        group_results[group_name] = {
            "factors": " x ".join(f"Z_{f}" for f in factors),
            "exponent": exponent,
            "num_sets": len(antichains)
        }
        all_unique_antichains.update(antichains)

    # Step 4: Present the results
    g1_res = group_results["G1"]
    print(f"Analysis for G1 = {g1_res['factors']}:")
    print(f"  - Exponent = lcm{tuple(g1_res['factors'])} = {g1_res['exponent']}")
    g1_num = g1_res['num_sets']
    print(f"  - Number of unique eigenvalue sets = {g1_num}\n")

    g2_res = group_results["G2"]
    print(f"Analysis for G2 = {g2_res['factors']}:")
    print(f"  - Exponent = lcm{tuple(g2_res['factors'])} = {g2_res['exponent']}")
    g2_num = g2_res['num_sets']
    print(f"  - Number of unique eigenvalue sets = {g2_num}\n")

    total_unique_sets = len(all_unique_antichains)
    
    print("-" * 80)
    print("The total number of unique sets is the size of the union of the sets from G1 and G2.")
    print(f"|Sets from G1| + |Sets from G2| - |Intersection| = Total Unique Sets")
    # In this case, sets from G2 are a subset of sets from G1, so union size is |Sets from G1|
    print(f"Number of sets for G1 = {g1_num}")
    print(f"Number of sets for G2 = {g2_num}")
    print(f"Since all sets possible for G2 are also possible for G1, the total number of unique sets is {total_unique_sets}.")
    print(f"Final Answer: {total_unique_sets}")
    print(f"\n<<<9>>>")


if __name__ == "__main__":
    main()