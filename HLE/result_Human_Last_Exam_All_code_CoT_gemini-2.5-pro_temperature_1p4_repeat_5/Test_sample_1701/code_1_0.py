from fractions import Fraction

def solve_bloom_xor_prob():
    """
    Calculates the probability that a bit in the derived Bloom filter B' is 1.
    """
    # Step 1: Define Bloom filter parameters
    m = 4  # Total bits in the Bloom filter B
    k = 2  # Number of hash functions
    n = 2  # Number of inputs

    kn = k * n

    print("Problem Parameters:")
    print(f"Total bits in B (m): {m}")
    print(f"Number of hash functions (k): {k}")
    print(f"Number of inputs (n): {n}")
    print(f"Total hash computations (kn): {kn}")
    print("-" * 30)

    # Step 2: Calculate the required probabilities using the Fraction class for precision
    # Probability that two specific bits (j and k) are both 0
    # P(b_j=0, b_k=0) = (1 - 2/m)^(kn)
    p_j0_k0_num = m - 2
    p_j0_k0_den = m
    p_j0_k0 = Fraction(p_j0_k0_num, p_j0_k0_den) ** kn

    # Probability that a single specific bit (j) is 0
    # P(b_j=0) = (1 - 1/m)^(kn)
    p_j0_num = m - 1
    p_j0_den = m
    p_j0 = Fraction(p_j0_num, p_j0_den) ** kn

    print("Intermediate Probability Calculations:")
    print(f"P(b_j=0) = (1 - 1/{m})^{kn} = ({p_j0_num}/{p_j0_den})^{kn} = {p_j0.numerator}/{p_j0.denominator}")
    print(f"P(b_j=0, b_k=0) = (1 - 2/{m})^{kn} = ({p_j0_k0_num}/{p_j0_k0_den})^{kn} = {p_j0_k0.numerator}/{p_j0_k0.denominator}")
    print("-" * 30)

    # Step 3: Calculate the joint probability P(b_j=1, b_k=0)
    # P(b_j=1, b_k=0) = P(b_k=0) - P(b_j=0, b_k=0)
    # By symmetry, P(b_k=0) is the same as P(b_j=0)
    p_j1_k0 = p_j0 - p_j0_k0

    # By symmetry, P(b_j=0, b_k=1) is the same
    p_j0_k1 = p_j1_k0

    # Step 4: Calculate the final probability
    # P(B'[i]=1) = P(b_j XOR b_k = 1) = P(b_j=1, b_k=0) + P(b_j=0, b_k=1)
    final_prob = p_j1_k0 + p_j0_k1

    print("Final Probability Calculation:")
    print(f"The probability of a bit in B' being 1 is given by P(b_j=1, b_k=0) + P(b_j=0, b_k=1).")
    print(f"P(b_j=1, b_k=0) = P(b_k=0) - P(b_j=0, b_k=0) = {p_j0.numerator}/{p_j0.denominator} - {p_j0_k0.numerator}/{p_j0_k0.denominator} = {p_j1_k0.numerator}/{p_j1_k0.denominator}")
    print(f"The final equation is: {p_j1_k0.numerator}/{p_j1_k0.denominator} + {p_j0_k1.numerator}/{p_j0_k1.denominator} = {final_prob.numerator}/{final_prob.denominator}")
    
    print("\nResult:")
    print(f"The probability that a randomly chosen bit in B' is 1 is {final_prob.numerator}/{final_prob.denominator}.")
    print(f"As a decimal, this is {float(final_prob)}.")


solve_bloom_xor_prob()