import sys

# In case n is large, we might need a higher recursion limit.
# For n=4, this is not strictly necessary.
sys.setrecursionlimit(2000)

# Memoization dictionaries to store computed values
memo_p = {}
memo_s = {}

def p(n):
    """
    Calculates the number of partitions of an integer n using a
    recursive approach with memoization.
    p_internal(n, k) calculates partitions of n using parts of size at most k.
    p(n) is equivalent to p_internal(n, n).
    """
    def p_internal(num, max_part):
        if num == 0:
            return 1
        if num < 0 or max_part <= 0:
            return 0
        if (num, max_part) in memo_p:
            return memo_p[(num, max_part)]
        
        # A partition of 'num' either uses 'max_part' as a part or it doesn't.
        # 1. Partitions not using 'max_part': p_internal(num, max_part - 1)
        # 2. Partitions using 'max_part': p_internal(num - max_part, max_part)
        result = p_internal(num, max_part - 1) + p_internal(num - max_part, max_part)
        memo_p[(num, max_part)] = result
        return result
        
    if n < 0: return 0
    return p_internal(n, n)


def S(n, k):
    """
    Calculates the Stirling number of the second kind, S(n, k),
    using the standard recurrence relation with memoization.
    S(n,k) is the number of ways to partition a set of n items into k non-empty subsets.
    """
    if (n, k) in memo_s:
        return memo_s[(n, k)]
    if k < 0 or k > n:
        return 0
    if n == 0 and k == 0:
        return 1
    if (n > 0 and k == 0) or (n == 0 and k > 0):
        return 0
    if k == 1 or k == n:
        return 1
    
    # Recurrence relation: S(n, k) = S(n-1, k-1) + k * S(n-1, k)
    result = S(n-1, k-1) + k * S(n-1, k)
    memo_s[(n, k)] = result
    return result

def solve():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size N.
    """
    N = 4
    
    print(f"The number of equivalence classes of endomorphisms on a set of size {N} is calculated by the formula:")
    print(f"Total = Sum_{{j=1 to {N}}} [ p(j) * S({N}, j) ]")
    print("\nwhere p(j) is the partition number and S(n, j) is the Stirling number of the second kind.")
    print("\nStep-by-step calculation:")

    total_classes = 0
    calculation_terms = []
    
    for j in range(1, N + 1):
        p_j = p(j)
        S_Nj = S(N, j)
        term = p_j * S_Nj
        total_classes += term
        
        print(f"\nFor j = {j}:")
        print(f"  - The number of partitions of {j}, p({j}), is {p_j}.")
        print(f"  - The Stirling number S({N}, {j}) is {S_Nj}.")
        print(f"  - The term for j={j} is p({j}) * S({N}, {j}) = {p_j} * {S_Nj} = {term}.")
        calculation_terms.append(str(term))
        
    print("\nFinal Summation:")
    final_equation = " + ".join(calculation_terms)
    print(f"Total = {final_equation} = {total_classes}")

    print(f"\nThus, there are {total_classes} distinct equivalence classes of endomorphisms on a set of size 4.")
    
solve()

<<<38>>>