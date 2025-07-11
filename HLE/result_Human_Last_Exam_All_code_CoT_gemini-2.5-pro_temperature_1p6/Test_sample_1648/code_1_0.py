import math

def solve_k_theory_problem(p, k):
    """
    This function solves the problem based on Mitchell's theorem on the vanishing of K-groups.
    
    The theorem states that the p-Sylow subgroup of K_{2n}(Z/p^k) is zero for
    n >= p^k / (p-1).
    
    The bound is sharp, so the largest n for which the group is non-zero is
    the largest integer less than this bound.
    """
    
    # The ring is Z/(p^k)
    ring_modulus = p**k
    
    # Calculate the vanishing bound for n
    bound = ring_modulus / (p - 1)
    
    # The largest integer n for which the group can be non-zero is
    # ceil(bound) - 1. For example, if bound is 13.5, n must be < 13.5.
    # The largest integer n is 13. ceil(13.5) is 14. 14 - 1 = 13.
    largest_n = math.ceil(bound) - 1
    
    print(f"The ring is Z/{ring_modulus}.")
    print(f"The parameters are p = {p}, k = {k}.")
    print(f"The vanishing condition for n is: n >= {p**k} / ({p}-1) = {bound}")
    print(f"Therefore, the p-Sylow subgroup of K_(2n)(Z/{ring_modulus}) is zero for all n >= {math.ceil(bound)}.")
    print(f"The largest natural number n for which the group is non-zero is {largest_n}.")

# Parameters for the specific problem Z/27
p = 3
k = 3

solve_k_theory_problem(p, k)