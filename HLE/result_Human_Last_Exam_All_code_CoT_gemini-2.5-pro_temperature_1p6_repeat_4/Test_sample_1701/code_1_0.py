import math
from fractions import Fraction

def solve_bloom_filter_probability():
    """
    Calculates the probability that a bit in the derived array B' is 1.
    """
    # Parameters of the Bloom filter
    m = 4  # number of bits in B
    k = 2  # number of hash functions
    n = 2  # number of inputs

    print(f"Bloom Filter Parameters: m={m}, k={k}, n={n}\n")

    # Step 1: Calculate the probability that a bit in B is 0 (p0)
    # p0 = (1 - 1/m)^(k*n)
    p0 = Fraction(m - 1, m) ** (k * n)
    
    print("The equation for the probability of a bit in B being 0 is: P(bit=0) = (1 - 1/m)^(k*n)")
    print(f"P(bit=0) = (1 - 1/{m})^({k}*{n}) = ({m-1}/{m})^{k*n} = {p0.numerator}/{p0.denominator}")
    print("-" * 20)

    # Step 2: Calculate the probability that a bit in B is 1 (p1)
    # p1 = 1 - p0
    p1 = 1 - p0
    
    print("The equation for the probability of a bit in B being 1 is: P(bit=1) = 1 - P(bit=0)")
    print(f"P(bit=1) = 1 - {p0.numerator}/{p0.denominator} = {p1.numerator}/{p1.denominator}")
    print("-" * 20)

    # Step 3: Calculate the probability that a bit in B' is 1
    # B'[i] = 1 if the two corresponding bits from B have different values.
    # P(B'[i]=1) = P(bit_a=0, bit_b=1) + P(bit_a=1, bit_b=0)
    # Since the bit states are independent, this is p0*p1 + p1*p0 = 2*p0*p1
    final_prob = 2 * p0 * p1

    print("The bit B'[i] is 1 if the corresponding bits from the two halves of B are different.")
    print("The equation for this probability is: P(B'[i]=1) = 2 * P(bit=0) * P(bit=1)")
    
    # Print the final equation with all numbers
    print(f"P(B'[i]=1) = 2 * {p0.numerator}/{p0.denominator} * {p1.numerator}/{p1.denominator}")
    print(f"             = {2 * p0.numerator * p1.numerator}/{p0.denominator * p1.denominator}")
    print(f"             = {final_prob.numerator}/{final_prob.denominator}")
    print("-" * 20)
    
    # Final answer
    print(f"\nThe final probability as a fraction is: {final_prob.numerator}/{final_prob.denominator}")
    print(f"As a decimal, the probability is approximately: {float(final_prob):.6f}")

solve_bloom_filter_probability()