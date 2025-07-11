import math

def solve_k_group_problem():
    """
    This function calculates the largest natural number n such that the (2n)th K-group of Z/27 is nonzero.
    The ring is Z/m where m = 27.
    This can be written as Z/p^k where p is a prime and k is a natural number.
    """
    m = 27

    # Find p and k such that m = p^k
    p = 0
    k = 0
    
    # We find the prime factor of m
    d = 2
    temp_m = m
    while d * d <= temp_m:
        if temp_m % d == 0:
            p = d
            break
        d += 1
    if p == 0:
        p = m

    # Now we find k
    if p != 0:
        k = int(round(math.log(m, p)))

    if p**k != m:
        print(f"The number {m} is not a perfect power of a prime number.")
        return

    print(f"The ring is Z/m where m = {m}. This is Z/p^k with p = {p} and k = {k}.")

    # A critical index in the study of K-groups of Z/p^k is given by N = (p^k - 1) / (p - 1).
    # It is known that K_{2n}(Z/p^k) vanishes for n >= N.
    # Therefore, the largest n for which the group is nonzero is N - 1.
    
    numerator = p**k - 1
    denominator = p - 1
    
    N = numerator // denominator
    
    largest_n = N - 1

    print(f"\nThe critical index N is calculated using the formula: N = (p^k - 1) / (p - 1)")
    print(f"N = ({p}^{k} - 1) / ({p} - 1)")
    print(f"N = ({p**k} - 1) / {denominator}")
    print(f"N = {numerator} / {denominator}")
    print(f"N = {N}")
    print(f"\nThe K-group is non-zero for n < N.")
    print(f"The largest natural number n is N - 1.")
    print(f"Largest n = {N} - 1 = {largest_n}")
    
    # Final answer
    print(f"\nFinal Answer: {largest_n}")

solve_k_group_problem()