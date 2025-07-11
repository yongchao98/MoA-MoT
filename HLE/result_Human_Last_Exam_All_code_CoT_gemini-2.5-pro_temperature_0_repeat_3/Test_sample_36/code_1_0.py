import math
from itertools import combinations

def get_divisors(n):
    """Returns the set of divisors of a positive integer n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def is_antichain(subset):
    """Checks if a subset of numbers is an antichain under divisibility."""
    # Singletons are always antichains.
    if len(subset) <= 1:
        return True
    # Check every pair for divisibility.
    for i in range(len(subset)):
        for j in range(i + 1, len(subset)):
            a, b = subset[i], subset[j]
            if a % b == 0 or b % a == 0:
                return False
    return True

def solve():
    """
    Calculates the number of unique sets S(rho) by counting the number of
    non-empty antichains in the poset of divisors of 18.
    """
    n = 18
    divisors = get_divisors(n)
    
    print(f"The Abelian groups of order 18 are C_18 and C_2 x C_3 x C_3.")
    print(f"The character orders for C_18 are the divisors of 18: {divisors}")
    print(f"The character orders for C_2 x C_3 x C_3 are the divisors of its exponent, 6: {get_divisors(6)}")
    print("The sets of eigenvalues from the second group are a subset of those from the first.")
    print("We count the number of unique sets by finding the number of non-empty antichains in the poset of divisors of 18.\n")

    total_antichains = 0
    max_len = len(divisors)
    
    # Iterate through all possible sizes of non-empty subsets
    for r in range(1, max_len + 1):
        count_for_size_r = 0
        # Generate all subsets of size r
        subsets_r = combinations(divisors, r)
        for subset in subsets_r:
            if is_antichain(subset):
                count_for_size_r += 1
        
        if count_for_size_r > 0:
            print(f"Number of antichains of size {r}: {count_for_size_r}")
        total_antichains += count_for_size_r

    print(f"\nTotal number of unique sets is the sum of these counts.")
    print(f"Final Answer: {total_antichains}")

solve()
<<<9>>>