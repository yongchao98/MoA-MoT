import fractions

def solve_bloom_xor_prob():
    """
    Calculates the probability that a bit in a derived Bloom filter is 1.

    The derived Bloom filter B' is created by XORing corresponding bits
    from the two halves of an original Bloom filter B.
    P(B'[i]=1) = P(B1[i] != B2[i])
                 = P(B1[i]=0, B2[i]=1) + P(B1[i]=1, B2[i]=0)
                 = 2 * P(bit in B is 0) * P(bit in B is 1)
    """
    # Parameters of the original Bloom filter B
    m = 4  # number of bits in B
    k = 2  # number of hash functions
    n = 2  # number of inputs

    print(f"Bloom Filter B has m={m} bits, k={k} hash functions, and n={n} inputs.")
    print("-" * 20)

    # 1. Calculate the probability that a bit in B is 0 (p0)
    # p0 = (1 - 1/m)^(k*n)
    p0_frac = fractions.Fraction(1 - 1, m) ** (k * n)
    p0_num = p0_frac.numerator
    p0_den = p0_frac.denominator

    print("Step 1: Calculate the probability (p0) that a bit in B is 0.")
    print(f"p0 = (1 - 1/{m}) ^ ({k} * {n})")
    print(f"p0 = (3/4) ^ 4")
    print(f"p0 = {p0_num}/{p0_den}")
    print()

    # 2. Calculate the probability that a bit in B is 1 (p1)
    # p1 = 1 - p0
    p1_frac = 1 - p0_frac
    p1_num = p1_frac.numerator
    p1_den = p1_frac.denominator

    print("Step 2: Calculate the probability (p1) that a bit in B is 1.")
    print(f"p1 = 1 - p0")
    print(f"p1 = 1 - {p0_num}/{p0_den}")
    print(f"p1 = {p1_num}/{p1_den}")
    print()

    # 3. Calculate the final probability for a bit in B' to be 1
    # P(B'[i]=1) = 2 * p0 * p1
    final_prob_frac = 2 * p0_frac * p1_frac
    final_prob_num = final_prob_frac.numerator
    final_prob_den = final_prob_frac.denominator

    print("Step 3: Calculate the probability P(B'[i]=1).")
    print("This is the probability that the corresponding bits in the two halves of B are different.")
    print("P(B'[i]=1) = 2 * p0 * p1")
    print(f"P(B'[i]=1) = 2 * ({p0_num}/{p0_den}) * ({p1_num}/{p1_den})")
    print(f"P(B'[i]=1) = {final_prob_num}/{final_prob_den}")
    print()
    
    # Print the final result as a decimal
    final_prob_decimal = float(final_prob_frac)
    print("Final answer as a decimal:")
    print(final_prob_decimal)


solve_bloom_xor_prob()
<<<14175/32768>>>