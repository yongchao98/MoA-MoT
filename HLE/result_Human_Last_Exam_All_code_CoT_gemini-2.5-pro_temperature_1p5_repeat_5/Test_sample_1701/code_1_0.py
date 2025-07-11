import fractions

def solve_bloom_filter_xor_prob():
    """
    Calculates the probability that a bit in a derived Bloom filter B' is 1.
    """
    # --- Parameters ---
    m = 4  # Total bits in the original Bloom filter B
    k = 2  # Number of hash functions
    n = 2  # Number of inputs

    # Total number of hash computations
    total_hashes = k * n

    # --- Calculations using fractions for precision ---
    # Convert numbers to Fraction objects
    frac_m = fractions.Fraction(m)
    frac_1 = fractions.Fraction(1)
    frac_2 = fractions.Fraction(2)

    # Probability that a specific bit B[j] remains 0
    # P(B[j]=0) = (1 - 1/m)^(k*n)
    p_b_is_0_term = (frac_1 - frac_1 / frac_m)
    p_b_is_0 = p_b_is_0_term ** total_hashes

    # Probability that two specific bits B[j] and B[l] both remain 0
    # P(B[j]=0, B[l]=0) = (1 - 2/m)^(k*n)
    p_b_both_0_term = (frac_1 - frac_2 / frac_m)
    p_b_both_0 = p_b_both_0_term ** total_hashes

    # The probability we want is P(B[j] != B[l]) = 2 * (P(B[j]=0) - P(B[j]=0, B[l]=0))
    p_diff = p_b_is_0 - p_b_both_0
    final_prob = 2 * p_diff

    # --- Print the step-by-step derivation ---
    print("This script calculates the probability P(B'[i]=1) where B'[i] = B[j] \u2295 B[l].")
    print("This is equivalent to P(B[j] != B[l]) = 2 * [P(B[j]=0) - P(B[j]=0, B[l]=0)].\n")
    print(f"Given parameters: m={m}, k={k}, n={n}")
    print(f"Total hash computations = k * n = {total_hashes}\n")

    print("The final calculation is:")
    print(f"P(B'[i]=1) = 2 * [ (1 - 1/{m})^({total_hashes}) - (1 - 2/{m})^({total_hashes}) ]")
    print(f"P(B'[i]=1) = 2 * [ ({p_b_is_0_term})^4 - ({p_b_both_0_term})^4 ]")
    print(f"P(B'[i]=1) = 2 * [ {p_b_is_0} - {p_b_both_0} ]")
    print(f"P(B'[i]=1) = 2 * [ {p_diff} ]")
    print(f"P(B'[i]=1) = {final_prob}\n")

    print(f"The final probability is {final_prob}.")
    print(f"As a decimal, this is {float(final_prob)}.")

solve_bloom_filter_xor_prob()