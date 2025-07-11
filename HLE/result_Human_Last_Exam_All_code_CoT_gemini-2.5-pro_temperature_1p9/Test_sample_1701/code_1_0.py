from fractions import Fraction

def solve_bloom_filter_probability():
    """
    Calculates the probability that a randomly chosen bit in a derived
    Bloom filter B' is equal to 1.
    """
    # Step 1: Define Bloom filter parameters
    m = 4  # number of bits in B
    k = 2  # number of hash functions
    n = 2  # number of inputs

    nk = n * k  # total number of hash computations

    print(f"Bloom Filter Parameters:")
    print(f"  - Size of B (m): {m}")
    print(f"  - Hash functions (k): {k}")
    print(f"  - Inserted items (n): {n}")
    print(f"  - Total hash lookups (n * k): {nk}\n")

    # Step 2: Define the goal
    print("Goal: Calculate P(B'[i] = 1)")
    print("This is equivalent to finding P(B[j] != B[j+2]) for j=1,2.")
    print("Let's calculate P(B'[1] = 1), which is P(B[1] != B[3]).\n")

    # Step 3: Formulate the probability condition
    print("Condition: B'[1] = 1 <=> (B[1]=1 and B[3]=0) or (B[1]=0 and B[3]=1).")
    print("P(B'[1]=1) = P(B[1]=1, B[3]=0) + P(B[1]=0, B[3]=1)\n")
    
    # Step 4: Calculate intermediate probabilities using the Fraction class for precision

    # P(B[j]=0) = probability that a specific bit is 0.
    # This occurs if all nk hashes miss this bit.
    # Probability of one hash missing is (m-1)/m.
    p_bit_is_0 = Fraction(m - 1, m) ** nk
    print(f"Calculating P(B[j]=0): ((m-1)/m)^(n*k) = (({m-1}/{m}))^{nk} = {p_bit_is_0}")

    # P(B[i]=0, B[j]=0) = probability that two specific bits are both 0.
    # This occurs if all nk hashes miss both bits.
    # Probability of one hash missing both is (m-2)/m.
    p_two_bits_are_0 = Fraction(m - 2, m) ** nk
    print(f"Calculating P(B[i]=0, B[j]=0): ((m-2)/m)^(n*k) = (({m-2}/{m}))^{nk} = {p_two_bits_are_0}\n")
    
    # P(B[1]=1, B[3]=0) = P(B[3]=0) - P(B[1]=0, B[3]=0)
    # This formula works because the event {B[3]=0} is the union of two disjoint events:
    # {B[1]=1, B[3]=0} and {B[1]=0, B[3]=0}.
    p_one_is_1_three_is_0 = p_bit_is_0 - p_two_bits_are_0
    print("Calculating the probability of one bit being 1 and another being 0:")
    print(f"P(B[1]=1, B[3]=0) = P(B[3]=0) - P(B[1]=0, B[3]=0)")
    print(f"P(B[1]=1, B[3]=0) = {p_bit_is_0} - {p_two_bits_are_0} = {p_one_is_1_three_is_0}\n")

    # By symmetry, P(B[1]=0, B[3]=1) is the same.
    p_one_is_0_three_is_1 = p_one_is_1_three_is_0
    print(f"By symmetry, P(B[1]=0, B[3]=1) is also {p_one_is_0_three_is_1}.\n")
    
    # Step 5: Final Calculation
    final_prob = p_one_is_1_three_is_0 + p_one_is_0_three_is_1
    print("The final probability is the sum of these two disjoint events:")
    print(f"P(B'[1]=1) = P(B[1]=1, B[3]=0) + P(B[1]=0, B[3]=1)")
    print(f"P(B'[1]=1) = {p_one_is_1_three_is_0} + {p_one_is_0_three_is_1} = {final_prob}")

    print(f"\nThe simplified final probability is: {final_prob.limit_denominator()}")
    print(f"As a decimal, this is: {float(final_prob)}")
    
solve_bloom_filter_probability()
<<<65/128>>>