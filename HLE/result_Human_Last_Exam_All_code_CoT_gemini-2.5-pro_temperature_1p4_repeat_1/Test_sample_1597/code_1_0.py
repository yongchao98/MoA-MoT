import collections

def get_prime_factor_set(n):
    """
    Returns the set of prime factors of n.
    Using frozenset so it can be a dictionary key.
    """
    factors = set()
    d = 2
    temp = n
    while d * d <= temp:
        if temp % d == 0:
            factors.add(d)
            while temp % d == 0:
                temp //= d
        d += 1
    if temp > 1:
        factors.add(temp)
    return frozenset(factors)

def solve():
    """
    Counts the number of allowed ordered pairs (a,b) with 1 <= a,b <= 1000.
    An ordered pair (a,b) is allowed if a=1, b=1, or for a,b > 1, P(a) = P(b).
    """
    limit = 1000

    # Case 1: a = 1 or b = 1
    count_case1 = limit + limit - 1

    # Case 2: a, b > 1 and P(a) = P(b)
    prime_factor_sets_counts = collections.defaultdict(int)
    for n in range(2, limit + 1):
        p_set = get_prime_factor_set(n)
        prime_factor_sets_counts[p_set] += 1

    count_case2 = 0
    for s in prime_factor_sets_counts:
        count = prime_factor_sets_counts[s]
        count_case2 += count * count

    total_count = count_case1 + count_case2
    
    # Print the final equation as requested.
    print(f"{count_case1} + {count_case2} = {total_count}")

solve()