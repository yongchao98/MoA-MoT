import fractions

def solve_probability():
    """
    Calculates and explains the probability that a bit in B' is 1.
    """
    # Problem parameters
    m = 4  # number of bits in B
    k = 2  # number of hash functions
    n = 2  # number of inputs

    # Total number of hash computations
    total_hashes = n * k

    # --- Calculations ---
    # To maintain precision, we will work with fractions.

    # 1. Calculate P(B[j]=0)
    # A single bit j is 0 if all 'total_hashes' miss this position.
    # The probability a single hash misses position j is (m-1)/m.
    # P(B[j]=0) = ((m-1)/m) ^ total_hashes
    p_bit_zero_num = (m - 1) ** total_hashes
    p_bit_zero_den = m ** total_hashes
    p_bit_zero = fractions.Fraction(p_bit_zero_num, p_bit_zero_den)

    # 2. Calculate P(B[j]=0, B[l]=0)
    # Two bits j and l are 0 if all 'total_hashes' miss both positions.
    # The probability a single hash misses both is (m-2)/m.
    # P(B[j]=0, B[l]=0) = ((m-2)/m) ^ total_hashes
    p_two_bits_zero_num = (m - 2) ** total_hashes
    p_two_bits_zero_den = m ** total_hashes
    p_two_bits_zero = fractions.Fraction(p_two_bits_zero_num, p_two_bits_zero_den)

    # 3. Calculate P(B[j]=0, B[l]=1) = P(B[j]=0) - P(B[j]=0, B[l]=0)
    p_j0_l1 = p_bit_zero - p_two_bits_zero

    # 4. Calculate the final probability P(B'[i]=1) = 2 * P(B[j]=0, B[l]=1)
    final_prob = 2 * p_j0_l1

    # --- Print the derivation with numbers ---
    print("The probability P(B'[i]=1) is derived as follows:")
    print("P(B'[i]=1) = P(B_1[i] != B_2[i])")
    print("           = P(B[j]=0, B[l]=1) + P(B[j]=1, B[l]=0)")
    print("           = 2 * (P(B[j]=0) - P(B[j]=0, B[l]=0))")
    print("\nWhere:")
    print(f"m = {m} (bits), n = {n} (inputs), k = {k} (hashes)")
    print(f"Total hashes = n * k = {n} * {k} = {total_hashes}")
    print("\nCalculations:")
    print(f"P(B[j]=0) = (({m}-1)/{m})^({n}*{k}) = ({m-1}/{m})^{total_hashes} = {p_bit_zero_num}/{p_bit_zero_den}")
    print(f"P(B[j]=0, B[l]=0) = (({m}-2)/{m})^({n}*{k}) = ({m-2}/{m})^{total_hashes} = {p_two_bits_zero_num}/{p_two_bits_zero_den}")
    print("\nFinal Equation:")
    print(f"P(B'[i]=1) = 2 * ({p_bit_zero_num}/{p_bit_zero_den} - {p_two_bits_zero_num}/{p_two_bits_zero_den})")
    print(f"           = 2 * (({p_bit_zero} - {p_two_bits_zero}))")
    print(f"           = 2 * ({p_j0_l1})")
    print(f"           = {final_prob}")
    print(f"\nThe final probability as a decimal is: {float(final_prob)}")

if __name__ == "__main__":
    solve_probability()