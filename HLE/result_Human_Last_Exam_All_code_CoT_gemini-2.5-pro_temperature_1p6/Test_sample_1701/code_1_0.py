import fractions

def solve_bloom_filter_probability():
    """
    Calculates the probability that a bit in a derived Bloom filter array is 1.
    """
    # --- Problem Parameters ---
    m = 4  # number of bits in the Bloom filter B
    k = 2  # number of hash functions
    n = 2  # number of inputs

    # Total number of hashes performed
    N = n * k

    print("Step 1: Understand the Goal")
    print("We want to find the probability P(B'[i]=1).")
    print(f"B' is defined such that B'[i] = B[j] XOR B[l] for some j, l.")
    print("This means B'[i]=1 if B[j] and B[l] have different values.")
    print("Due to symmetry, we only need to calculate this for one pair of bits, e.g., B[1] and B[3].")
    print("P(B'[i]=1) = P(B[j] != B[l]) = P(B[j]=0, B[l]=1) + P(B[j]=1, B[l]=0)")
    print("-" * 30)

    print("Step 2: Probability of a single bit being 0 in B")
    print("A bit B[j] is 0 if none of the N = n*k hash operations land on its position.")
    p_bit_is_zero_num = (m - 1)**N
    p_bit_is_zero_den = m**N
    p_bit_is_zero = fractions.Fraction(p_bit_is_zero_num, p_bit_is_zero_den)

    print(f"Total hashes N = n * k = {n} * {k} = {N}")
    print(f"P(B[j]=0) = (({m}-1)/{m})^{N} = {p_bit_is_zero.numerator}/{p_bit_is_zero.denominator}")
    print("-" * 30)

    print("Step 3: Probability of two specific bits being 0 in B")
    print("Bits B[j] and B[l] (j!=l) are both 0 if none of the N hashes land on either position.")
    p_two_bits_are_zero_num = (m - 2)**N
    p_two_bits_are_zero_den = m**N
    p_two_bits_are_zero = fractions.Fraction(p_two_bits_are_zero_num, p_two_bits_are_zero_den)

    print(f"P(B[j]=0, B[l]=0) = (({m}-2)/{m})^{N} = {p_two_bits_are_zero.numerator}/{p_two_bits_are_zero.denominator}")
    print("-" * 30)

    print("Step 4: Calculate the final probability")
    print("We need P(B[j]=0, B[l]=1). We know P(B[j]=0) = P(B[j]=0, B[l]=0) + P(B[j]=0, B[l]=1).")
    print("Therefore, P(B[j]=0, B[l]=1) = P(B[j]=0) - P(B[j]=0, B[l]=0).")
    p_01 = p_bit_is_zero - p_two_bits_are_zero

    print(f"P(B[j]=0, B[l]=1) = {p_bit_is_zero} - {p_two_bits_are_zero} = {p_01}")
    print("\nBy symmetry, P(B[j]=1, B[l]=0) is the same.")
    
    print("\nFinally, P(B'[i]=1) = P(B[j]=0, B[l]=1) + P(B[j]=1, B[l]=0).")
    final_prob = p_01 + p_01
    
    # --- Final Answer Calculation ---
    print("\n--- Summary of Final Equation ---")
    print(f"The probability that a randomly chosen bit in B' is 1 is:")
    print(f"P(B'[i]=1) = 2 * ( P(B[j]=0) - P(B[j]=0, B[l]=0) )")
    print(f"           = 2 * ( (({m}-1)/{m})^{N} - (({m}-2)/{m})^{N} )")
    
    common_den = p_bit_is_zero_den
    p_two_bits_common_den_num = p_two_bits_are_zero.numerator * (common_den // p_two_bits_are_zero.denominator)
    p_01_common_den = fractions.Fraction(p_bit_is_zero.numerator - p_two_bits_common_den_num, common_den)
    
    print(f"           = 2 * ( {p_bit_is_zero.numerator}/{common_den} - {p_two_bits_common_den_num}/{common_den} )")
    print(f"           = 2 * ( {p_01_common_den.numerator}/{p_01_common_den.denominator} )")
    print(f"           = {final_prob.numerator}/{final_prob.denominator}")

    print(f"\nThe probability as a fraction is {final_prob}.")
    print(f"The probability as a decimal is {float(final_prob)}.")


solve_bloom_filter_probability()