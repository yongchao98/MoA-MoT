import math

def calculate_probability():
    """
    Calculates the probability that for each individual, there exists a type of item
    for which they hold strictly more copies than any other individual.
    """

    # Step 1: Define combinatorial helper numbers
    # D(n) is the number of derangements of n items.
    def D(n):
        if n < 0: return 0
        if n == 0: return 1
        return n * D(n-1) + (-1)**n

    D5 = D(5)  # 44
    num_5_cycles = math.factorial(5 - 1) // 2  # 12
    num_3_2_cycles_config = math.comb(5, 3)  # 10

    # Step 2: Calculate N_fav_fixed, the number of favorable sequences for a fixed specialty assignment.
    # We sum the number of sequences for each case based on the diagonal value d.
    
    # Case d=5: C = 5*I
    # Each person's hand is (5,0,0,0,0). Arrangements per hand = 5!/5! = 1
    # Number of such matrices = 1
    sequences_d5 = 1 * (1**5)

    # Case d=4: C = 4*I + P_p, where p is a derangement
    # Each hand is a permutation of (4,1,0,0,0). Arrangements per hand = 5!/4! = 5
    # Number of such matrices = D5
    sequences_d4 = D5 * (5**5)

    # Case d=3:
    # Subcase 3a: Off-diagonal is 2s (derangement). Hand is (3,2,0,0,0).
    # Arrangements per hand = 5!/(3!2!) = 10. Number of matrices = D5.
    sequences_d3a = D5 * (10**5)
    
    # Subcase 3b: Off-diagonal is 1s (5-cycle). Hand is (3,1,1,0,0).
    # Arrangements per hand = 5!/(3!1!1!) = 20. Number of matrices = 12.
    sequences_d3b = num_5_cycles * (20**5)
    
    # Subcase 3c: Mixed 1s and 2s (3-cycle + pair).
    # Two hands are (3,2,0,0,0) type (arr=10), three are (3,1,1,0,0) type (arr=20).
    # Arrangements per matrix = 10^2 * 20^3. Number of matrices = 10.
    sequences_d3c = num_3_2_cycles_config * (10**2 * 20**3)

    N_fav_fixed = sequences_d5 + sequences_d4 + sequences_d3a + sequences_d3b + sequences_d3c

    # Step 3: Calculate total favorable arrangements F
    # Multiply by 5! for all possible specialty-to-person assignments.
    F = math.factorial(5) * N_fav_fixed

    # Step 4: Calculate total arrangements S
    S = math.factorial(25) // (math.factorial(5)**5)

    # Step 5: Calculate probability P = F/S and simplify the fraction
    common_divisor = math.gcd(F, S)
    numerator = F // common_divisor
    denominator = S // common_divisor

    print(f"Total number of ways to distribute the items (S): {S}")
    print(f"Number of favorable distributions (F): {F}")
    print(f"The probability is P = F/S")
    print(f"P = {F} / {S}")
    print(f"P = {numerator} / {denominator}")


calculate_probability()