def solve(n, k):
    """
    Calculates the k-th element of S_n using the optimal bitwise method.
    The value is S_n[k] = v(k+1) + 1, where v(x) is the number of trailing zeros in x.
    This implementation prints the components of the calculation as requested.
    """
    # According to the discovered pattern, the value is independent of n.
    # We only need k to compute the result.
    x = k + 1

    # This algorithm is valid for k >= 0, so x >= 1.
    if x <= 0:
        # This case should not be reached given the problem constraints.
        print("Error: k must be non-negative.")
        return

    # Efficiently find the number of trailing zeros, v(x).
    # 1. Isolate the least significant bit (LSB). The result is a power of 2.
    lsb = x & -x
    # 2. The number of trailing zeros is the base-2 logarithm of the LSB.
    v_x = lsb.bit_length() - 1

    # The value of the element is v(x) + 1.
    result = v_x + 1

    # As requested, output the numbers in the final equation: v(k+1) + 1 = result
    print(f"{v_x} + 1 = {result}")


# Example from the problem description: n = 2, k = 3. Expected output: 3
# Calculation: k+1=4. v(4)=2. Result = 2+1=3. The print should be "2 + 1 = 3"
solve(2, 3)