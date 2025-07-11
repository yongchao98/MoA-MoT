import fractions

def solve_bloom_filter_xor_probability():
    """
    Calculates the probability that a bit in a derived XOR'd Bloom filter array is 1.
    """
    # Problem parameters
    m = 4  # Total bits in the original Bloom filter B
    k = 2  # Number of hash functions
    n = 2  # Number of inputs

    # Total number of hash calculations
    N = k * n

    print("### Problem Setup ###")
    print(f"Original Bloom filter B has m={m} bits.")
    print(f"We have n={n} inputs and k={k} hash functions.")
    print(f"Total hash calculations N = k * n = {N}.")
    print("B' is derived such that B'[i] = B[i] XOR B[i+m/2].")
    print("We want to find P(B'[i] = 1), which is P(B[i] != B[i+m/2]).\n")
    
    # We choose i=1 and j=3 for our analysis without loss of generality.
    # P(B[1] != B[3]) = 1 - P(B[1] = B[3])
    # P(B[1] = B[3]) = P(B[1]=0, B[3]=0) + P(B[1]=1, B[3]=1)

    print("### Calculation Steps ###")

    # Step 1: Calculate P(B[1]=0 and B[3]=0)
    # This occurs if none of the N hashes map to 1 or 3.
    # Prob for one hash to miss both is (m-2)/m.
    p_b0_and_c0_num = (m - 2)**N
    p_b0_and_c0_den = m**N
    p_b0_and_c0 = fractions.Fraction(p_b0_and_c0_num, p_b0_and_c0_den)
    print(f"1. P(B[1]=0 and B[3]=0) = (({m}-2)/{m})^{N} = ({m-2}/{m})^{N} = {p_b0_and_c0}")

    # Step 2: Calculate P(B[1]=1 and B[3]=1)
    # We use the principle of inclusion-exclusion: P(A and C) = 1 - (P(not A) + P(not C) - P(not A and not C))
    # P(not A) is P(B[1]=0). Prob for one hash to miss 1 is (m-1)/m.
    p_b0_num = (m - 1)**N
    p_b0_den = m**N
    p_b0 = fractions.Fraction(p_b0_num, p_b0_den)
    p_c0 = p_b0 # by symmetry

    # P(not A or not C) = P(B[1]=0) + P(B[3]=0) - P(B[1]=0, B[3]=0)
    p_b0_or_c0 = p_b0 + p_c0 - p_b0_and_c0

    # P(B[1]=1 and B[3]=1) = 1 - P(not A or not C)
    p_b1_and_c1 = 1 - p_b0_or_c0
    print(f"2. To find P(B[1]=1 and B[3]=1), we first find the probability of its complement, P(B[1]=0 or B[3]=0).")
    print(f"   P(B[1]=0) = (({m}-1)/{m})^{N} = ({m-1}/{m})^{N} = {p_b0}")
    print(f"   P(B[1]=0 or B[3]=0) = P(B[1]=0) + P(B[3]=0) - P(B[1]=0 and B[3]=0)")
    print(f"   = {p_b0} + {p_c0} - {p_b0_and_c0} = {p_b0_or_c0}")
    print(f"   P(B[1]=1 and B[3]=1) = 1 - {p_b0_or_c0} = {p_b1_and_c1}")

    # Step 3: Calculate P(B[1] = B[3])
    p_b_equals_c = p_b0_and_c0 + p_b1_and_c1
    print(f"3. P(B[1]=B[3]) = P(B[1]=0, B[3]=0) + P(B[1]=1, B[3]=1)")
    print(f"   = {p_b0_and_c0} + {p_b1_and_c1} = {p_b_equals_c}")

    # Step 4: Calculate the final probability P(B[1] != B[3])
    final_prob = 1 - p_b_equals_c
    print(f"4. The final probability is 1 - P(B[1]=B[3])")
    print(f"   = 1 - {p_b_equals_c} = {final_prob}\n")

    print("### Final Answer ###")
    print(f"The probability that a randomly chosen bit in B' is 1 is {final_prob}.")
    print(f"As a fraction, this is {final_prob.numerator}/{final_prob.denominator}.")

solve_bloom_filter_xor_probability()