import math

def solve_proportionality_puzzle():
    """
    This script solves the proportionality puzzle by calculating the minimum
    number of voters (s1, s2) required to leave voter 1 unsatisfied
    under PJR and EJR, respectively.
    """

    # --- Problem Parameters ---
    committee_size_k = 100
    
    # We consider the critical group G = {1, 2, 3, 4, 5, 6}, which is present
    # in any extended profile.
    critical_group_size = 6
    
    # The intersection of approved candidates for this group G is {a,b,c}.
    intersection_size = 3
    
    print("--- Calculating s1 for Proportional Justified Representation (PJR) ---")
    
    # For voter 1 to be unsatisfied, their approved candidates {a,b,c,x} cannot be in the committee.
    # This means the group G's common intersection {a,b,c} cannot be in the committee.
    # PJR guarantees representation for G if |G| >= n/k.
    # To avoid this guarantee, we must have |G| < n/k.
    
    print(f"The critical group has size |G| = {critical_group_size}.")
    print(f"The committee size is k = {committee_size_k}.")
    print("PJR guarantees representation if the following inequality holds:")
    print(f"|G| >= n / k   =>   {critical_group_size} >= n / {committee_size_k}")
    
    # This simplifies to n <= k * |G|
    pjr_bound = committee_size_k * critical_group_size
    print(f"This is equivalent to: n <= {committee_size_k} * {critical_group_size}, which means n <= {pjr_bound}.")
    print("\nTo allow voter 1 to be unsatisfied, this guarantee must be avoided.")
    print(f"Therefore, we must have n > {pjr_bound}.")
    
    # s1 is the smallest integer n such that n > 600.
    s1 = pjr_bound + 1
    print(f"The smallest integer profile size s1 is {s1}.\n")
    
    
    print("--- Calculating s2 for Extended Justified Representation (EJR) ---")

    # EJR guarantees representation for G if |G| >= l * n/k for any l <= |intersection|.
    # The intersection size is 3, so we check for l = 1, 2, 3.
    # To avoid the guarantee, we must have |G| < l * n/k for ALL applicable l.
    # This is equivalent to n > k * |G| / l for all applicable l.
    
    print(f"The critical group has size |G| = {critical_group_size} and intersection size {intersection_size}.")
    print("EJR provides a guarantee if |G| >= l * n/k for l in {1, 2, 3}.")

    bounds_for_n = []
    for l in range(1, intersection_size + 1):
        # We need n > k * |G| / l
        ejr_bound_for_l = math.floor(committee_size_k * critical_group_size / l)
        bounds_for_n.append(ejr_bound_for_l)
        print(f"\nFor l = {l}:")
        print(f"EJR is triggered if: {critical_group_size} >= {l} * n / {committee_size_k}")
        print(f"This simplifies to n <= ({committee_size_k} * {critical_group_size}) / {l}, which is n <= {ejr_bound_for_l}.")
        print(f"To avoid this, we need n > {ejr_bound_for_l}.")
        
    # To avoid the EJR guarantee completely for this group, n must be larger than the
    # maximum of all these bounds.
    n_must_be_greater_than = max(bounds_for_n)
    print("\nTo avoid the EJR guarantee for all relevant l, n must satisfy the strictest condition.")
    print(f"The strictest condition is n > {n_must_be_greater_than}.")

    # s2 is the smallest integer n such that n > 600.
    s2 = n_must_be_greater_than + 1
    print(f"The smallest integer profile size s2 is {s2}.\n")

    print("--- Final Result ---")
    final_pair = (s1, s2)
    print(f"The solution pair (s1, s2) is: {final_pair}")

solve_proportionality_puzzle()