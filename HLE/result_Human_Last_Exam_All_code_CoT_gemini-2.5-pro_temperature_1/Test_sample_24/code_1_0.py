import math

def n_subgroups(i, j, n):
    """Calculates one term in the summation part of the formula."""
    # N(C2, i) and N(C5, j) are 1 for the relevant i and j, so we omit them.
    # The formula for the term is: n * (i+j-1)! / (i! * j!)
    numerator = n * math.factorial(i + j - 1)
    denominator = math.factorial(i) * math.factorial(j)
    return numerator // denominator

def solve():
    """
    Calculates the number of subgroups of index 7 in C2 * C5.
    """
    n = 7
    
    # N(C2, 7) is 0 since 7 does not divide 2.
    N_A_n = 0
    # N(C5, 7) is 0 since 7 does not divide 5.
    N_B_n = 0
    
    print(f"The number of subgroups of index {n} in C2 is {N_A_n}.")
    print(f"The number of subgroups of index {n} in C5 is {N_B_n}.")
    
    # The pairs (i, j) for which N(C2, i) and N(C5, j) are non-zero are:
    # (1, 1), (1, 5), (2, 1), (2, 5)
    
    # i=1, j=1
    term1 = n_subgroups(1, 1, n)
    print(f"Contribution for (i=1, j=1): {n} * ({1}+{1}-1)! / ({1}! * {1}!) = {term1}")
    
    # i=1, j=5
    term2 = n_subgroups(1, 5, n)
    print(f"Contribution for (i=1, j=5): {n} * ({1}+{5}-1)! / ({1}! * {5}!) = {term2}")
    
    # i=2, j=1
    term3 = n_subgroups(2, 1, n)
    print(f"Contribution for (i=2, j=1): {n} * ({2}+{1}-1)! / ({2}! * {1}!) = {term3}")

    # i=2, j=5
    term4 = n_subgroups(2, 5, n)
    print(f"Contribution for (i=2, j=5): {n} * ({2}+{5}-1)! / ({2}! * {5}!) = {term4}")

    total_subgroups = N_A_n + N_B_n + term1 + term2 + term3 + term4
    
    print("\nFinal Calculation:")
    print(f"Total = {N_A_n} + {N_B_n} + {term1} + {term2} + {term3} + {term4} = {total_subgroups}")
    
    print(f"\nThe total number of subgroups of index 7 in G = C2 * C5 is {total_subgroups}.")

solve()