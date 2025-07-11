import math

def solve_sequence(n, k):
    """
    Calculates the k-th element of the sequence S_n using an optimal method.

    The sequence S_n is defined as:
    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1}

    The value of S_n[k] can be determined directly by the formula:
    S_n[k] = v_2(k + 1) + 1
    where v_2(x) is the 2-adic valuation of x, which is the exponent of 2 in
    the prime factorization of x, or equivalently, the number of trailing
    zeros in the binary representation of x.

    This is the most optimal implementation as its complexity depends only on the
    number of bits in k, not on the value of n.
    """
    if k < 0:
        print("Error: Index k cannot be negative.")
        return

    # Use 1-based indexing for the calculation, as the pattern is clearer.
    p = k + 1

    # An efficient way to find the number of trailing zeros is to isolate the
    # least significant bit (LSB) and find its log base 2.
    # p & -p isolates the LSB.
    # The number of trailing zeros is log2 of the LSB.
    # In Python, for a power of two `x`, math.log2(x) is equivalent to
    # x.bit_length() - 1.
    if p == 0:
        # This case happens if k = -1, handled by the check above.
        # If it were possible, the value would be undefined.
        num_trailing_zeros = -1 # As a convention
    else:
        lsb = p & -p
        num_trailing_zeros = lsb.bit_length() - 1
    
    result = num_trailing_zeros + 1

    # The prompt requires outputting each number in the final equation.
    # The final equation is the relationship between the inputs and the output.
    print(f"S_{n}[{k}] = {result}")

# Example from the problem description
# n=2, k=3 should result in 3
# solve_sequence(2, 3)
# S_2 = [1, 2, 1, 3, 1, 2, 1]
# S_2[3] = 3

# Another example: n=2, k=5
# S_2[5] = 2
# solve_sequence(2, 5)

# Example with larger numbers, fulfilling the constraints n, k > INT_MAX/2
# For example, INT_MAX is 2147483647
n_large = 2147483648
k_large = 3458764513820540927 # A 62-bit number
solve_sequence(n_large, k_large)
