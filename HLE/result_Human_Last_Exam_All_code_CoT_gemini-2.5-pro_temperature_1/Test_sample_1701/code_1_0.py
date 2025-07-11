import math

def solve_bloom_filter_xor_prob():
    """
    Calculates the probability that a bit in the derived Bloom filter B' is 1.
    """
    # Parameters of the Bloom filter B
    m = 4  # number of bits in B
    k = 2  # number of hash functions
    n = 2  # number of inputs

    # Total number of hash values generated
    nk = n * k
    
    # The common denominator for our fractions will be m^nk
    common_denominator = m**nk

    # Step 1: Calculate the probability that a single specific bit in B is 0.
    # This occurs if none of the nk hashes land on this bit.
    # For a single hash, the probability of not hitting the bit is (m-1)/m.
    # For nk hashes, this is ((m-1)/m)^nk.
    prob_b_j_is_0_num = (m - 1)**nk
    
    print(f"The setup involves a Bloom filter with m={m} bits, k={k} hash functions, and n={n} inputs.")
    print(f"The probability of a specific bit B[j] being 0 is P(B[j]=0) = (({m}-1)/{m})^({n*k}) = {prob_b_j_is_0_num}/{common_denominator}")

    # Step 2: Calculate the probability that two specific, distinct bits in B are both 0.
    # This occurs if none of the nk hashes land on either of these two bits.
    # For a single hash, the probability is (m-2)/m.
    # For nk hashes, this is ((m-2)/m)^nk.
    prob_b_j0_l0_num = (m - 2)**nk
    
    print(f"The probability of two specific bits B[j] and B[l] both being 0 is P(B[j]=0, B[l]=0) = (({m}-2)/{m})^({n*k}) = {prob_b_j0_l0_num}/{common_denominator}")

    # Step 3: Calculate the probability P(B'[i]=1).
    # P(B'[i]=1) = P(B[j] != B[l]) = 2 * P(B[j]=0, B[l]=1)
    # We know P(B[j]=0, B[l]=1) = P(B[j]=0) - P(B[j]=0, B[l]=0)
    prob_b_j0_l1_num = prob_b_j_is_0_num - prob_b_j0_l0_num
    
    # The final probability is 2 * P(B[j]=0, B[l]=1)
    final_prob_num = 2 * prob_b_j0_l1_num
    final_prob_den = common_denominator

    # Simplify the fraction
    common_divisor = math.gcd(final_prob_num, final_prob_den)
    simplified_num = final_prob_num // common_divisor
    simplified_den = final_prob_den // common_divisor

    print("\nThe probability we want to find is P(B'[i]=1), which equals 2 * (P(B[j]=0) - P(B[j]=0, B[l]=0)).")
    print(f"P(B'[i]=1) = 2 * ({prob_b_j_is_0_num}/{common_denominator} - {prob_b_j0_l0_num}/{common_denominator})")
    print(f"P(B'[i]=1) = 2 * ({prob_b_j0_l1_num}/{common_denominator})")
    print(f"P(B'[i]=1) = {final_prob_num}/{final_prob_den}")
    print(f"\nThe final probability is {final_prob_num}/{final_prob_den}, which simplifies to {simplified_num}/{simplified_den} or as a decimal {final_prob_num/final_prob_den}.")

solve_bloom_filter_xor_prob()