import math

def get_prime_factors(n):
    """
    Calculates the set of distinct prime factors of a given integer n.
    """
    factors = set()
    if n == 1:
        return frozenset()
    
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
    Counts the number of allowed pairs (a,b) with 1 <= a,b <= 1000.
    An ordered pair (a,b) is allowed if P(a) is a subset of P(b) or P(b) is a subset of P(a).
    """
    limit = 1000
    
    # Step 1 & 2: Group numbers by their set of prime factors
    prime_sets_map = {}
    for i in range(1, limit + 1):
        s = get_prime_factors(i)
        prime_sets_map.setdefault(s, 0)
        prime_sets_map[s] += 1
        
    # Convert map to a list of (set, count) items for easier iteration
    sets_with_counts = list(prime_sets_map.items())
    
    num_allowed_pairs = 0
    
    # Step 3: Iterate through all pairs of groups and check the condition
    for i in range(len(sets_with_counts)):
        for j in range(len(sets_with_counts)):
            s_i, n_i = sets_with_counts[i]
            s_j, n_j = sets_with_counts[j]
            
            # Check if one set of prime factors is a subset of the other
            if s_i.issubset(s_j) or s_j.issubset(s_i):
                # If so, all pairs of numbers from these two groups are allowed
                num_allowed_pairs += n_i * n_j
                
    print(num_allowed_pairs)

solve()