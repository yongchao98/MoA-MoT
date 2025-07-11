import fractions

def solve_bloom_filter_probability():
    """
    Calculates the probability that a bit in a derived bit array B' is 1.

    The problem is defined by a Bloom filter B with 4 bits, 2 hash functions,
    and 2 inputs. B is split into B1 and B2. B' is the XOR of B1 and B2.
    We want to find P(B'[i] = 1).
    """
    # Parameters from the problem description
    m = 4  # Number of bits in the Bloom filter
    k = 2  # Number of hash functions
    n = 2  # Number of inputs inserted

    # Total number of hash operations
    nk = n * k

    # The probability we want to find is P(B'[i] = 1). Due to symmetry, we can
    # calculate this for i=0: P(B'[0] = 1) = P(B[0] XOR B[2] = 1).
    # This is equal to P(B[0] != B[2]).
    # The final probability, P, can be derived as P = 2 * (A - B), where:
    # A = P(B[j]=0), the probability a single bit is 0.
    # B = P(B[j]=0, B[l]=0), the probability two distinct bits are both 0.

    print("This problem can be solved by calculating the probabilities of bit states in the Bloom filter.")
    print("The final probability P is given by the equation: P = 2 * (A - B), where A and B are defined as follows:\n")

    # --- Calculation of intermediate probabilities A and B ---
    
    # Calculate A = P(B[j]=0)
    # A bit is 0 if all nk hashes miss it. A single hash misses with probability (m-1)/m.
    prob_A_num = (m - 1)**nk
    prob_A_den = m**nk
    prob_A = fractions.Fraction(prob_A_num, prob_A_den)

    # Calculate B = P(B[j]=0, B[l]=0)
    # Two bits are 0 if all nk hashes miss both. A single hash misses both with probability (m-2)/m.
    prob_B_num = (m - 2)**nk
    prob_B_den = m**nk
    prob_B = fractions.Fraction(prob_B_num, prob_B_den)

    # --- Outputting the numbers for the final equation ---
    print("--- Numbers for the Final Equation ---")
    
    # Print the value for A
    print(f"A = P(a single bit is 0) = (({m}-1)/{m})^({nk}) = {prob_A.numerator}/{prob_A.denominator}")

    # Print the value for B
    print(f"B = P(two distinct bits are 0) = (({m}-2)/{m})^({nk}) = {prob_B.numerator}/{prob_B.denominator}\n")

    # --- Final Calculation and Output ---

    # To show the subtraction clearly, represent B with the same denominator as A
    prob_B_common_den = prob_B.limit_denominator(prob_A.denominator)

    # The difference (A - B)
    diff_prob = prob_A - prob_B
    
    # The final probability P = 2 * (A - B)
    final_prob = 2 * diff_prob

    print("--- The Final Equation Solved Step-by-Step ---")
    print(f"P = 2 * (A - B)")
    print(f"P = 2 * ({prob_A.numerator}/{prob_A.denominator} - {prob_B_common_den.numerator}/{prob_B_common_den.denominator})")
    print(f"P = 2 * ({diff_prob.numerator}/{diff_prob.denominator})")
    print(f"P = {final_prob.numerator}/{final_prob.denominator}\n")
    
    print("--- Final Answer ---")
    print(f"The probability that for a randomly chosen bit in B', B'[i]=1, is {final_prob}.")
    print(f"As a decimal, this is equal to {float(final_prob)}.")

if __name__ == '__main__':
    solve_bloom_filter_probability()