def solve():
    """
    Calculates the size of the largest union of 20 antichains in the divisor poset of N = 2^2 * 3^30.
    This is equivalent to counting the number of divisors x of N such that Omega(x) <= 19.
    For x = 2^a * 3^b, Omega(x) = a+b.
    The constraints are: 0 <= a <= 2, 0 <= b <= 30, and a+b <= 19.
    """
    max_a = 2
    max_omega = 19
    
    total_count = 0
    terms = []
    
    # Iterate through possible exponents of 2, which are a = 0, 1, 2
    for a in range(max_a + 1):
        # For a given 'a', we need a+b <= 19, so b <= 19-a.
        # Since b must also be non-negative (b >= 0),
        # b can range from 0 to 19-a.
        # The number of possible integer values for b is (19-a) - 0 + 1.
        
        # We also have the constraint b <= 30, but since 19-a is always <= 19, this constraint is always met.
        
        max_b = max_omega - a
        
        # Number of values for b is from 0 to max_b, inclusive.
        count_for_a = max_b + 1
        
        total_count += count_for_a
        terms.append(str(count_for_a))

    equation = " + ".join(terms)
    print(f"The size of the largest union of 20 antichains is found by counting the number of divisors 2^a * 3^b of N = 2^2 * 3^30 that satisfy a+b <= 19.")
    print(f"The number of such divisors is the sum of counts for each possible value of a (from 0 to 2).")
    print(f"The final calculation is: {equation} = {total_count}")

solve()