import math

def solve():
    """
    Solves the problem of decomposing 10^100.
    
    The logic is as follows:
    1.  The problem describes a unique way to represent any number as a sum of terms
        of the form tet(2, a) * pow(2, b). This uniqueness implies the number of sequences is 1.
        Each term tet(2, a) * pow(2, b) is equal to 2^(2^(a-1) + b).
        The constraint tet(2, a) > pow(2, b) implies b < 2^(a-1).
        This structure means that any power of two, 2^k, has a unique representation
        in this form, derived from its binary expansion.
    
    2.  The "largest multiplication in the sequence", T = tet(2, a1) * pow(2, b1), corresponds
        to the largest power of 2 that is part of the sum for N = 10^100.
        This largest power of 2 is 2^k_max, where k_max = floor(log2(10^100)).
    
    3.  We calculate k_max = floor(100 * log2(10)).
    
    4.  From k_max, we determine a1 and b1. The term T is 2^k_max, which must be equal to
        2^(2^(a1-1) + b1). So, k_max = 2^(a1-1) + b1.
    
    5.  Given the constraint b1 < 2^(a1-1), a1 is found by satisfying 2^(a1-1) <= k_max < 2^a1.
        This gives a1 - 1 = floor(log2(k_max)).
    
    6.  b1 is then calculated as the remainder: b1 = k_max - 2^(a1-1).
    """
    
    # Based on the analysis, the number of sequences is 1.
    count = 1
    
    # Calculate k_max for the largest term
    # k_max = floor(log2(10^100)) = floor(100 * log2(10))
    k_max = math.floor(100 * math.log2(10))
    
    # Calculate a1 for k_max
    # We need 2^(a1-1) <= k_max < 2^a1
    # This means a1-1 = floor(log2(k_max))
    a1_minus_1 = math.floor(math.log2(k_max))
    a1 = a1_minus_1 + 1
    
    # Calculate b1 for k_max
    # b1 = k_max - 2^(a1-1)
    b1 = k_max - (2**a1_minus_1)
    
    print(f"{count} {a1} {b1}")

solve()