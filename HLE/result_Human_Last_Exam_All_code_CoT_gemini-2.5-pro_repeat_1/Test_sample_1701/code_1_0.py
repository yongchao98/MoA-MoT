import fractions

def solve_bloom_filter_probability():
    """
    Calculates the probability that a bit in a derived Bloom filter array is 1.

    The problem defines a Bloom filter B with m=4 bits, k=2 hash functions,
    and n=2 inputs. The total number of hash computations is N = n * k = 4.

    The filter B is split into B1=B[1,2] and B2=B[3,4]. A new array B' is
    derived such that B'[i] = B1[i] XOR B2[i]. We want to find the probability
    that a randomly chosen bit in B' is 1.

    This is equivalent to finding P(B[1] != B[3]) or P(B[2] != B[4]).
    Due to symmetry, these probabilities are equal.

    The calculation follows these steps:
    1. P(B[i]!=B[j]) = P(B[i]=1, B[j]=0) + P(B[i]=0, B[j]=1)
    2. By symmetry, this is 2 * P(B[i]=1, B[j]=0).
    3. P(B[i]=1, B[j]=0) = P(B[j]=0) - P(B[i]=0, B[j]=0).
    4. Calculate P(B[j]=0) and P(B[i]=0, B[j]=0) based on the N hash computations.
    """
    # Problem parameters
    m = 4  # bits in the Bloom filter B
    k = 2  # hash functions
    n = 2  # inputs

    # Total number of hash computations
    N = n * k

    print("Problem Parameters:")
    print(f"  Number of bits in B (m): {m}")
    print(f"  Number of hash functions (k): {k}")
    print(f"  Number of inputs (n): {n}")
    print(f"  Total hash computations (N = n * k): {N}")
    print("-" * 40)

    print("Step-by-step Calculation:")
    print("Let's find P(B'[i]=1) = P(B_1[i] XOR B_2[i] = 1).")
    print("This is equivalent to finding P(B[x] != B[y]) for a corresponding pair of bits.")
    print("-" * 40)

    # Step 1: Calculate P(B[j]=0)
    # The probability that a single hash does NOT land on a specific bit is (m-1)/m.
    # The probability that all N hashes do NOT land on that bit is ((m-1)/m)^N.
    p_j_is_0_num = m - 1
    p_j_is_0_den = m
    prob_j_is_0 = fractions.Fraction(p_j_is_0_num, p_j_is_0_den) ** N
    print(f"1. Probability that a specific bit is 0, e.g., P(B[3]=0):")
    print(f"   P(B[3]=0) = ( (m-1) / m ) ^ N")
    print(f"   P(B[3]=0) = ( ({p_j_is_0_num}) / {p_j_is_0_den} ) ^ {N}")
    print(f"   P(B[3]=0) = {prob_j_is_0.numerator}/{prob_j_is_0.denominator}")
    print()

    # Step 2: Calculate P(B[i]=0, B[j]=0)
    # The probability that a single hash misses two specific bits is (m-2)/m.
    # The probability that all N hashes miss both bits is ((m-2)/m)^N.
    p_ij_is_0_num = m - 2
    p_ij_is_0_den = m
    prob_ij_is_0 = fractions.Fraction(p_ij_is_0_num, p_ij_is_0_den) ** N
    print(f"2. Probability that two specific bits are 0, e.g., P(B[1]=0, B[3]=0):")
    print(f"   P(B[1]=0, B[3]=0) = ( (m-2) / m ) ^ N")
    print(f"   P(B[1]=0, B[3]=0) = ( ({p_ij_is_0_num}) / {p_ij_is_0_den} ) ^ {N}")
    print(f"   P(B[1]=0, B[3]=0) = {prob_ij_is_0.numerator}/{prob_ij_is_0.denominator}")
    print()

    # Step 3: Calculate P(B[i]=1, B[j]=0)
    # P(B[i]=1, B[j]=0) = P(B[j]=0) - P(B[i]=0, B[j]=0)
    prob_i1_j0 = prob_j_is_0 - prob_ij_is_0
    print(f"3. Probability that one bit is 1 and the other is 0, e.g., P(B[1]=1, B[3]=0):")
    print(f"   P(B[1]=1, B[3]=0) = P(B[3]=0) - P(B[1]=0, B[3]=0)")
    print(f"   P(B[1]=1, B[3]=0) = {prob_j_is_0} - {prob_ij_is_0}")
    print(f"   P(B[1]=1, B[3]=0) = {prob_i1_j0.numerator}/{prob_i1_j0.denominator}")
    print()

    # Step 4: Calculate the final probability P(B[i] != B[j])
    # P(B[i] != B[j]) = 2 * P(B[i]=1, B[j]=0)
    final_prob = 2 * prob_i1_j0
    print(f"4. Final probability that the XOR is 1, P(B[1] != B[3]):")
    print(f"   P(B'[i]=1) = P(B[1] != B[3]) = 2 * P(B[1]=1, B[3]=0)")
    print(f"   P(B'[i]=1) = 2 * {prob_i1_j0}")
    print(f"   P(B'[i]=1) = {final_prob.numerator}/{final_prob.denominator}")
    print("-" * 40)

    print(f"The final probability is {final_prob}.")
    print(f"As a decimal, this is {float(final_prob)}.")

solve_bloom_filter_probability()