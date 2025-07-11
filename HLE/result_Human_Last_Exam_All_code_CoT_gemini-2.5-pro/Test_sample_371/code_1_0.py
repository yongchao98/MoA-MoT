import sys

def solve():
    """
    Calculates the k-th element of the sequence S_n, where S_0 = [1]
    and S_n = S_{n-1} + [n+1] + S_{n-1}.

    This function implements the optimal solution for large n and k.
    """
    # The problem specifies n, k > INT_MAX/2. We'll use example values
    # that satisfy this constraint. INT_MAX is 2**31 - 1.
    # We choose n and k to be large integers.
    try:
        n = 2 * (10**9)
        k = 3 * (10**18) + 1234567890
    except OverflowError:
        print("Running on a 32-bit system, using smaller example numbers.", file=sys.stderr)
        n = 30
        k = 1000

    # The k-th element of S_n follows the formula: S_n[k] = v_2(k + 1) + 1
    # where v_2(x) is the number of trailing zeros in the binary representation of x.
    # The value of n is a constraint on the valid range of k (k < 2**(n+1) - 1),
    # but is not needed for the calculation itself.

    # We first calculate k + 1. Python's ints handle arbitrary size.
    x = k + 1

    # To find the number of trailing zeros of a positive integer x in Python,
    # we can use the expression (x & -x).bit_length() - 1.
    # The operation `x & -x` isolates the least significant bit (LSB) of x.
    # For a power of two, `p`, `p.bit_length()` returns `log2(p) + 1`.
    # The number of trailing zeros is `log2` of the LSB, hence the -1.
    if x <= 0:
        # This case should not be reached with the problem constraints (k>=0)
        trailing_zeros = -1
        result = 0
    else:
        trailing_zeros = (x & -x).bit_length() - 1
        result = trailing_zeros + 1

    # The prompt requests to output each number in the final equation.
    # The final calculation is `result = trailing_zeros + 1`.
    print(f"{trailing_zeros} + 1 = {result}")

solve()