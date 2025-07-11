import math

def count_divisors(n):
    """
    Calculates tau(n), the number of divisors of n.
    """
    if n == 0:
        return 0
    count = 0
    i = 1
    while i * i <= n:
        if n % i == 0:
            # If i is a divisor, n/i is also a divisor.
            if i * i == n:
                count += 1
            else:
                count += 2
        i += 1
    return count

def solve():
    """
    Calculates the number of allowed pairs (a,b) with 1 <= a,b <= 1000.
    The condition for a pair to be allowed is a|b or b|a.
    """
    limit = 1000

    # Calculate S = Sum_{k=1 to limit} tau(k)
    # This is equivalent to the number of pairs (a,b) such that a|b.
    s = 0
    for k in range(1, limit + 1):
        s += count_divisors(k)

    # By symmetry, the number of pairs with b|a is also S.
    # The number of pairs with a=b is 'limit'.
    # Using inclusion-exclusion principle:
    # Total = (a|b) + (b|a) - (a=b)
    # Total = S + S - limit
    total_allowed_pairs = 2 * s - limit
    
    print(f"The number of pairs (a,b) with a|b is {s}.")
    print(f"The number of pairs with b|a is also {s}.")
    print(f"The number of pairs with a=b is {limit}.")
    print(f"By the Principle of Inclusion-Exclusion, the total number of allowed pairs is 2 * {s} - {limit} = {total_allowed_pairs}.")

solve()
