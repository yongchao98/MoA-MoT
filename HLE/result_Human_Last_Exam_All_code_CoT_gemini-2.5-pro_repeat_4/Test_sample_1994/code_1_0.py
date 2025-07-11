import math

def solve_tetration_sum():
    """
    Calculates the number of sequences and the a1, b1 values for the largest term
    in the tetration-based sum representation of 10^100.
    """
    
    # The number of sequences is always 1 because the representation is unique.
    count = 1
    
    # The target number is N = 10^100.
    # The largest term in the sum corresponds to the most significant bit of N.
    # Its exponent, E_max, is floor(log2(N)).
    # E_max = floor(log2(10^100)) = floor(100 * log2(10))
    E_max = math.floor(100 * math.log2(10))
    
    # Now, we find the (a1, b1) pair corresponding to E_max.
    # The relationship is E_max = 2^(a1 - 1) + b1, with the constraint b1 < 2^(a1 - 1).
    # This implies 2^(a1 - 1) <= E_max < 2^a1.
    
    # From 2^(a1 - 1) <= E_max, we take log2: a1 - 1 <= log2(E_max).
    # So, a1 - 1 = floor(log2(E_max)).
    a1 = math.floor(math.log2(E_max)) + 1
    
    # Now we find b1 using the relationship E_max = 2^(a1 - 1) + b1.
    # We use integer exponentiation (1 << x is equivalent to 2**x for integers).
    power_of_2 = 1 << (a1 - 1)
    b1 = E_max - power_of_2
    
    # The problem asks for three numbers: the count of sequences, a1, and b1.
    print(f"{count} {a1} {b1}")

solve_tetration_sum()