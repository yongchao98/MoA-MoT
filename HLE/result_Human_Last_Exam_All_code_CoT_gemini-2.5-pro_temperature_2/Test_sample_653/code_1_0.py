import math
from itertools import permutations

def solve_snp_order():
    """
    Solves for the SNP order based on the provided genetic expression data.
    The function demonstrates that no solution exists with the given numbers under the standard model.
    """
    
    # The five aFC values provided in the problem.
    # aFC = Mutant Allele Expression / WT Allele Expression
    afc_values = [1/3, 1/2, 3/2, 2, 3]
    afc_ranks = {1/3: 1, 1/2: 2, 3/2: 3, 2: 4, 3: 5}

    # As derived in the plan, the problem implies a multiplicative model for gene expression.
    # The product of all aFCs is a key term in our equations.
    F_prod = math.prod(afc_values)

    # We derived two equations for the aFCs at the ends of the SNP block, x_end1 and x_end2.
    # Equation 1: From the individual with expression equal to the WT level.
    # F_prod + x_end1 = 2
    x_end1 = 2 - F_prod

    # Equation 2: From the individual with expression 2.5 times the WT level.
    # F_prod + x_end2 = 5
    x_end2 = 5 - F_prod
    
    print("Derived model equations:")
    print(f"Product of all aFCs (F_prod) = {F_prod:.2f}")
    print("Equation for first endpoint (x_end1): F_prod + x_end1 = 2.0")
    print(f"--> {F_prod:.2f} + {x_end1:.2f} = 2.0")
    print("Equation for second endpoint (x_end2): F_prod + x_end2 = 5.0")
    print(f"--> {F_prod:.2f} + {x_end2:.2f} = 5.0\n")

    print(f"Solving these equations gives required endpoint aFCs: {x_end1:.2f} and {x_end2:.2f}")

    # Check if these calculated values exist in the provided list of aFCs.
    x_end1_found = any(math.isclose(x_end1, val) for val in afc_values)
    x_end2_found = any(math.isclose(x_end2, val) for val in afc_values)

    print(f"Provided aFCs: {sorted(afc_values)}")
    print(f"Is {x_end1:.2f} in the list? {'Yes' if x_end1_found else 'No'}")
    print(f"Is {x_end2:.2f} in the list? {'Yes' if x_end2_found else 'No'}\n")

    if x_end1_found and x_end2_found:
        print("A solution may exist. But since one endpoint must be 3.5, no solution is found.")
    else:
        print("Conclusion: No combination of the given aFCs can satisfy the conditions of the problem.")
        print("This suggests an inconsistency in the problem statement (e.g., a typo in the aFC values or the expression level ratio).")
        # Let's assume 2.5 was a typo for 2.25, which gives x_end2=3.
        # This makes the endpoints 1/2 and 3. Ranks 2 and 5.
        # The middle ranks are 1, 3, 4. The problem provides no information to order them.
        # The problem is underdetermined. As no unique answer can be found, we cannot provide a specific order.
        # A possible solution under this assumption would be 21345.
        
solve_snp_order()
# A specific order is requested despite the inconsistency. Based on a plausible typo correction (2.5x -> 2.25x),
# the endpoints become aFC=0.5 (rank 2) and aFC=3 (rank 5).
# The internal order {1/3, 3/2, 2} (ranks 1,3,4) is not determined by the problem statement.
# To provide a single answer as requested, a plausible simple order for the internal part is ascending rank.
# This results in the order 21345.
# Let's provide this as the most reasonable guess under the circumstances.
print("\nSince a single answer is required, we assume a typo in the problem (2.5 -> 2.25).")
print("This makes the endpoints Rank 2 and Rank 5. Arranging the middle ranks (1,3,4) by value gives:")
print("21345")
#<<<21345>>>