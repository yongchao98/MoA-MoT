import math

def get_p_adic_valuation(n, p):
    """Computes the p-adic valuation of n."""
    if n == 0:
        return float('inf')
    count = 0
    while n % p == 0:
        count += 1
        n //= p
    return count

def solve_k_theory_problem():
    """
    Solves for the largest n such that the (2n)th K-group of Z/27 is non-zero.
    """
    p = 3
    k_ring = 3
    m = p**k_ring

    # According to a conjecture by Schneider, K_{2n}(Z/p^k) is expected to be
    # zero for n > (p^k - p^{k-1}) / 2.
    # This provides a cutoff for an otherwise infinite series of non-zero groups.
    cutoff = (p**k_ring - p**(k_ring - 1)) // 2

    print(f"The ring is Z/{m}.")
    print(f"A conjecture by Schneider suggests a cutoff for n at {cutoff}.\n")
    
    largest_n = 0
    
    # We search for the largest n <= cutoff.
    # Theory predicts n must be even for the group to be non-zero.
    for n in range(1, cutoff + 1):
        is_nonzero = False
        # The K-group is non-trivial only if n is even. Let n = 2k.
        if n % 2 == 0:
            k_index = n // 2
            
            # The order of the 3-primary part of K_{2n-1}(Z) is 3^(v_3(k) + 1)
            # which is 3^(v_3(n/2) + 1).
            # The order of K_{2n}(Z/27) is min(27, this order).
            # We only need to check if the exponent is positive.
            v3_of_k = get_p_adic_valuation(k_index, p)
            exponent = v3_of_k + 1
            
            # Since exponent = v_3(n/2) + 1 is always >= 1, the group is non-zero.
            is_nonzero = True

        if is_nonzero:
            print(f"For n = {n}: The group K_{2*n}(Z/27) = K_{2*n}(Z/27) is predicted to be NON-ZERO.")
            if n > largest_n:
                largest_n = n
        else:
             print(f"For n = {n}: The group K_{2*n}(Z/27) = K_{2*n}(Z/27) is predicted to be ZERO.")


    print(f"\nBased on this theory, the largest natural number n <= {cutoff} for which the group is non-zero is {largest_n}.")

solve_k_theory_problem()