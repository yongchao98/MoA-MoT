from fractions import Fraction

def solve_bloom_filter_xor_prob():
    """
    Calculates the probability that a bit in the derived Bloom filter B' is 1.
    
    The problem is defined by:
    - Bloom filter B with m=4 bits, k=2 hash functions, n=2 inputs.
    - B is split into B1 (bits 1,2) and B2 (bits 3,4).
    - A new filter B' is defined by B'[i] = B1[i] XOR B2[i].
    
    The probability P(B'[i]=1) is the probability that B1[i] and B2[i] are different.
    This is calculated as 1 - P(B'[i]=0), which is 1 - [P(B1[i]=0, B2[i]=0) + P(B1[i]=1, B2[i]=1)].
    """
    
    # Parameters
    m = 4  # size of bloom filter
    k = 2  # number of hash functions
    n = 2  # number of inputs
    kn = k * n

    # Let j and l be two distinct bit positions in the original Bloom filter B.
    # These correspond to B1[i] and B2[i] for some i.
    
    # Step 1: Calculate P(B[j]=0 and B[l]=0)
    # This happens if all kn hashes miss both bits j and l.
    # The probability for one hash to miss both is (m-2)/m.
    p_both_zero = Fraction(m - 2, m) ** kn

    # Step 2: Calculate P(B[j]=1 and B[l]=1) using inclusion-exclusion.
    # P(B[j]=1, B[l]=1) = 1 - P(B[j]=0 or B[l]=0)
    # P(B[j]=0 or B[l]=0) = P(B[j]=0) + P(B[l]=0) - P(B[j]=0, B[l]=0)
    
    # P(B[j]=0): probability that a specific bit is 0.
    # All kn hashes must miss bit j. Prob for one hash to miss is (m-1)/m.
    p_one_zero = Fraction(m - 1, m) ** kn
    
    # P(B[j]=0 or B[l]=0)
    p_zero_union = p_one_zero + p_one_zero - p_both_zero
    
    # P(B[j]=1 and B[l]=1)
    p_both_one = 1 - p_zero_union
    
    # Step 3: Calculate the probability that B'[i]=0 (the bits are the same)
    p_xor_is_zero = p_both_zero + p_both_one
    
    # Step 4: The final probability P(B'[i]=1) is 1 minus the probability that B'[i]=0
    final_prob = 1 - p_xor_is_zero
    
    # Print the explanation and the step-by-step equation
    print("The goal is to calculate P(B'[i] = 1), where B'[i] = B1[i] XOR B2[i].")
    print("This is equivalent to 1 - P(B'[i] = 0).")
    print("P(B'[i] = 0) occurs if B1[i] and B2[i] are the same (both 0 or both 1).\n")

    print(f"Let B[j] and B[l] be the two distinct bits in the original filter corresponding to B1[i] and B2[i].")
    print(f"P(B[j]=0, B[l]=0) = ((m-2)/m)^(kn) = (({m-2}/{m}))^{kn} = {p_both_zero}")
    print(f"P(B[j]=1, B[l]=1) = 1 - [P(B[j]=0) + P(B[l]=0) - P(B[j]=0, B[l]=0)]")
    print(f"                   = 1 - [({p_one_zero}) + ({p_one_zero}) - ({p_both_zero})] = {p_both_one}\n")
    
    print("The full calculation for the final probability is:")
    print(f"P(B'[i] = 1) = 1 - (P(B[j]=0, B[l]=0) + P(B[j]=1, B[l]=1))")
    print(f"             = 1 - ({p_both_zero} + {p_both_one})")
    print(f"             = 1 - {p_xor_is_zero}")
    print(f"             = {final_prob}\n")
    print(f"The final probability is {final_prob.numerator}/{final_prob.denominator}")

solve_bloom_filter_xor_prob()