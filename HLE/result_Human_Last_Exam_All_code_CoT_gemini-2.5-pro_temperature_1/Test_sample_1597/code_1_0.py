def solve():
    """
    Calculates the number of allowed ordered pairs (a, b) of positive integers
    where 1 <= a, b <= 1000.

    A pair (a, b) is allowed if for every primitive function f, f(ab) >= f(a)f(b).
    This condition is equivalent to the set of prime factors of a being a subset
    of the set of prime factors of b, or vice-versa.
    Let P(n) be the set of prime factors of n. The condition is P(a) ⊆ P(b) or P(b) ⊆ P(a).
    """
    MAX_N = 1000

    # Step 1: Generate sets of prime factors for all numbers up to MAX_N.
    # p_factors[i] will store the set of prime factors of i.
    # We use a sieve-like method for efficiency.
    p_factors = [set() for _ in range(MAX_N + 1)]
    for i in range(2, MAX_N + 1):
        if not p_factors[i]:  # This means i is a prime number
            for j in range(i, MAX_N + 1, i):
                p_factors[j].add(i)

    # Step 2: Count the pairs (a, b) that satisfy the condition.
    allowed_pairs_count = 0
    for a in range(1, MAX_N + 1):
        for b in range(1, MAX_N + 1):
            set_a = p_factors[a]
            set_b = p_factors[b]
            
            # Check if one set of prime factors is a subset of the other.
            if set_a.issubset(set_b) or set_b.issubset(set_a):
                allowed_pairs_count += 1
    
    # Step 3: Print the final count.
    print(allowed_pairs_count)

solve()