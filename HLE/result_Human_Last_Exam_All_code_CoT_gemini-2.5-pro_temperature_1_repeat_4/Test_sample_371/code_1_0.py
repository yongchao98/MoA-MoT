import sys

def find_kth_element_in_Sn(n, k):
    """
    Calculates the k-th element of the sequence S_n, where S_0 = [1] and
    S_n = S_{n-1} + [n+1] + S_{n-1}.

    The optimal solution is based on the observation that S_n[k] is equal to
    the number of trailing zeros in the binary representation of (k+1), plus one.
    This value is also known as v_2(k+1) + 1.

    Args:
        n (int): The sequence number.
        k (int): The index in the sequence (starting from 0).
    """

    # The problem states n and k can be very large, but the logic holds.
    # Python's arbitrary-precision integers handle this seamlessly.
    if k < 0:
        print("Error: Index k must be a non-negative integer.", file=sys.stderr)
        return

    # The value to inspect is k+1.
    x = k + 1

    # We calculate the number of trailing zeros in x's binary representation.
    # An efficient method is to isolate the lowest set bit and find its position.
    # `x & -x` in two's complement arithmetic isolates the lowest set bit.
    # The result is a power of 2, say 2^p, where p is the number of trailing zeros.
    # For a number y = 2^p, its bit length is p + 1.
    # Thus, p = y.bit_length() - 1.
    lowest_set_bit = x & -x
    num_trailing_zeros = lowest_set_bit.bit_length() - 1

    # The value of the element S_n[k] is num_trailing_zeros + 1.
    result = num_trailing_zeros + 1

    # Per the instructions, we output the numbers in the final equation.
    # The equation is: result = num_trailing_zeros + 1.
    # We present it in the context of the problem.
    print(f"S_{n}[{k}] = {num_trailing_zeros} + 1 = {result}")

if __name__ == '__main__':
    # Use the example from the problem description: n=2, k=3.
    # S_2 = [1, 2, 1, 3, 1, 2, 1], so S_2[3] should be 3.
    # Let's verify with our function.
    n_example = 2
    k_example = 3
    find_kth_element_in_Sn(n_example, k_example)