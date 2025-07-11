def solve():
    """
    Solves the problem for the two given primes.
    """

    # Define the primes
    p1 = 80039
    p2 = 80077

    # Case 1: p = 80039
    # The value is F(3) = 7/8 mod p1
    val1_num = 7
    val1_den = 8
    inv_den1 = pow(val1_den, -1, p1)
    ans1 = (val1_num * inv_den1) % p1
    
    # Case 2: p = 80077
    # The value is F(N) = 16 * F(7) = 16 * (27/64) = 27/4 mod p2
    val2_num = 27
    val2_den = 4
    inv_den2 = pow(val2_den, -1, p2)
    ans2 = (val2_num * inv_den2) % p2

    print(f"For p = {p1}, the expression to calculate is: ({val1_num} * pow({val1_den}, -1, {p1})) mod {p1}")
    print(f"For p = {p2}, the expression to calculate is: ({val2_num} * pow({val2_den}, -1, {p2})) mod {p2}")
    
    print(f"The results are: {ans1},{ans2}")

solve()