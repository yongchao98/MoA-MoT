from fractions import Fraction

def solve_bloom_filter_xor_prob():
    """
    Calculates the probability that a bit in a derived Bloom filter B' is 1.

    The problem is defined by a Bloom filter B with:
    - m = 4 bits
    - k = 2 hash functions
    - n = 2 inputs

    B is split into B1=B[1,2] and B2=B[3,4].
    B' is derived as B'[i] = B1[i] XOR B2[i].
    We want to find P(B'[i] = 1).
    """

    # Parameters
    m = 4  # total bits in B
    k = 2  # number of hash functions
    n = 2  # number of inputs

    print("This script calculates the probability that a randomly chosen bit in B' is 1.")
    print("The formula used is: P(B'[i]=1) = 2 * (P(a specific bit in B is 0) - P(two specific bits in B are 0))\n")
    print("Parameters:")
    print(f"  - Size of B (m): {m}")
    print(f"  - Number of hash functions (k): {k}")
    print(f"  - Number of inputs (n): {n}\n")

    # --- Step 1: Calculate P(B[j]=0) ---
    # This is the probability that a single specific bit in B is 0.
    # Formula: P(B[j]=0) = (1 - 1/m)^(n*k)
    hashes = n * k
    p_one_bit_zero_num = (m - 1)**hashes
    p_one_bit_zero_den = m**hashes
    p_one_bit_zero = Fraction(p_one_bit_zero_num, p_one_bit_zero_den)

    print("Step 1: Calculate the probability that a specific bit in B is 0, P(B[j]=0).")
    print(f"P(B[j]=0) = (1 - 1/{m})^({n}*{k}) = (({m-1})/{m})^{hashes} = {p_one_bit_zero_num}/{p_one_bit_zero_den}")
    print(f"P(B[j]=0) = {p_one_bit_zero}\n")

    # --- Step 2: Calculate P(B[j]=0, B[l]=0) ---
    # This is the probability that two specific bits in B (one from each half) are 0.
    # Formula: P(B[j]=0, B[l]=0) = (1 - 2/m)^(n*k)
    p_two_bits_zero_num = (m - 2)**hashes
    p_two_bits_zero_den = m**hashes
    p_two_bits_zero = Fraction(p_two_bits_zero_num, p_two_bits_zero_den)

    print("Step 2: Calculate the probability that two specific bits in B are 0, P(B[j]=0, B[l]=0).")
    print(f"P(B[j]=0, B[l]=0) = (1 - 2/{m})^({n}*{k}) = (({m-2})/{m})^{hashes} = {p_two_bits_zero_num}/{p_two_bits_zero_den}")
    print(f"P(B[j]=0, B[l]=0) = {p_two_bits_zero}\n")

    # --- Step 3: Calculate the final probability P(B'[i]=1) ---
    # Formula: 2 * (P(B[j]=0) - P(B[j]=0, B[l]=0))
    diff = p_one_bit_zero - p_two_bits_zero
    final_prob = 2 * diff

    print("Step 3: Calculate the final probability P(B'[i]=1).")
    print("The final equation is:")
    print(f"P(B'[i]=1) = 2 * (P(B[j]=0) - P(B[j]=0, B[l]=0))")
    print(f"P(B'[i]=1) = 2 * ({p_one_bit_zero} - {p_two_bits_zero})")
    print(f"P(B'[i]=1) = 2 * ({diff.numerator}/{diff.denominator})")
    print(f"P(B'[i]=1) = {final_prob.numerator}/{final_prob.denominator}\n")

    print("Final Answer:")
    print(f"The probability is {final_prob}, which is approximately {float(final_prob):.5f}.")

if __name__ == '__main__':
    solve_bloom_filter_xor_prob()
    # The final probability is 65/128
    # The submission format requires <<<answer>>>
    print("<<<65/128>>>")
