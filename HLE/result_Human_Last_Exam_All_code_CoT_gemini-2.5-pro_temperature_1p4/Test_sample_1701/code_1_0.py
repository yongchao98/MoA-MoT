import math

def solve_bloom_filter_xor_probability():
    """
    Calculates the probability that a randomly chosen bit in B' is 1.
    B' is derived from a Bloom filter B by XORing bits from its two halves.
    """
    # Parameters of the Bloom filter
    m = 4  # number of bits in the filter B
    k = 2  # number of hash functions
    n = 2  # number of inputs

    # Total number of hashing operations
    H = n * k

    print(f"Given a Bloom filter B with m={m} bits, k={k} hash functions, and n={n} inputs.")
    print(f"The total number of hash computations is H = n * k = {H}.")
    print(f"B is split into B1=B[1,2] and B2=B[3,4]. A new filter B' is formed by B'[i] = B1[i] XOR B2[i].")
    print(f"We want to find the probability that a random bit in B' is 1, i.e., P(B[j] XOR B[k] = 1) for j in B1 and k in B2.\n")

    # Step 1: Probability that a specific bit B[i] is 0
    # P(B[i]=0) = ((m-1)/m)^H
    prob_bit_is_0_num = (m - 1)**H
    prob_bit_is_0_den = m**H

    # Step 2: Probability that two specific bits B[i] and B[j] are both 0
    # P(B[i]=0, B[j]=0) = ((m-2)/m)^H
    prob_two_bits_are_0_num = (m - 2)**H
    prob_two_bits_are_0_den = m**H
    
    # Step 3: We need P(B[i]=0, B[j]=1), which is P(B[i]=0) - P(B[i]=0, B[j]=0)
    prob_i0_j1_num = prob_bit_is_0_num - prob_two_bits_are_0_num
    prob_i0_j1_den = prob_bit_is_0_den

    # Step 4: The final probability P(B[i] XOR B[j] = 1) is 2 * P(B[i]=0, B[j]=1)
    final_prob_num = 2 * prob_i0_j1_num
    final_prob_den = prob_i0_j1_den
    
    # Simplify the final fraction
    common_divisor = math.gcd(final_prob_num, final_prob_den)
    simple_num = final_prob_num // common_divisor
    simple_den = final_prob_den // common_divisor

    # Final result calculation string
    print("The probability P(B'[i] = 1) is given by the formula:")
    print("P = 2 * (P(A specific bit is 0) - P(Two specific bits are 0))")
    print("P = 2 * ( ( (m-1)/m )^H - ( (m-2)/m )^H )")
    print(f"P = 2 * ( ( ({m}-1)/{m} )^{H} - ( ({m}-2)/{m} )^{H} )")
    print(f"P = 2 * ( ({prob_bit_is_0_num}/{prob_bit_is_0_den}) - ({prob_two_bits_are_0_num}/{prob_two_bits_are_0_den}) )")
    print(f"P = 2 * ( {prob_i0_j1_num}/{prob_i0_j1_den} )")
    print(f"P = {final_prob_num}/{final_prob_den}")
    print(f"P = {simple_num}/{simple_den}")
    print(f"P = {final_prob_num / final_prob_den}")

solve_bloom_filter_xor_probability()
<<<65/128>>>