import sys

def solve_kth_element(n, k):
    """
    Calculates the k-th element of the sequence S_n, where:
    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1} for n > 0

    The optimal implementation is based on the observation that the value S_n[k]
    is equal to v_2(k+1) + 1, where v_2(x) is the 2-adic valuation of x,
    i.e., the number of trailing zeros in the binary representation of x.

    This function calculates this value using efficient bitwise operations.
    The value of n is not used in the calculation, as the result only depends on k.
    This is valid as long as k is a valid index in S_n (i.e., k < 2**(n+1) - 1),
    which is implied by the problem statement.
    """
    if not isinstance(n, int) or not isinstance(k, int) or k < 0:
        print("Error: n and k must be integers, and k must be non-negative.", file=sys.stderr)
        return

    # The value is based on k+1
    x = k + 1

    # We can find the number of trailing zeros by isolating the
    # least significant bit (LSB). In two's complement arithmetic,
    # the expression (x & -x) results in a number that is a power of 2,
    # representing the LSB of x.
    lsb = x & -x

    # For a number 'p' that is a power of two, its 0-indexed bit position
    # (which equals the number of trailing zeros) can be found with p.bit_length() - 1.
    num_trailing_zeros = lsb.bit_length() - 1

    # The final result is the number of trailing zeros plus one.
    result = num_trailing_zeros + 1

    # Per the instructions, output each number in the final equation.
    # The equation is: S_n[k] = v_2(k+1) + 1
    # which we calculate as: result = num_trailing_zeros + 1
    print(f"{result} = {num_trailing_zeros} + 1")


# Example from the problem description to demonstrate the function
n_example = 2
k_example = 3
solve_kth_element(n_example, k_example)
