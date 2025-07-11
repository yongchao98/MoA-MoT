import sys

def solve(n, k):
    """
    Calculates the k-th element of the sequence S_n.
    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1}

    The value S_n[k] is given by the formula v_2(k+1) + 1, where v_2(x) is
    the number of trailing zeros in the binary representation of x.
    """
    
    # The value of n does not affect the result, only the maximum valid k.
    # We check if k is a valid index for S_n.
    # The length of S_n is 2**(n+1) - 1.
    if k >= (1 << (n + 1)) - 1:
        print(f"Error: k={k} is out of bounds for S_{n} which has length {(1 << (n + 1)) - 1}.")
        return

    x = k + 1

    # To find the number of trailing zeros (p), we can isolate the
    # least significant bit (LSB) of x.
    # The expression (x & -x) in two's complement arithmetic results in a power of 2,
    # specifically 2**p, which has only the LSB of x set.
    # For example, if x = 12 (binary 1100), x & -x = 4 (binary 100).
    # The bit_length() of a power of 2, 2**p, is p+1.
    # So, p = (2**p).bit_length() - 1.
    
    if x == 0:
        # This case is not expected since k >= 0.
        p = -1 # Undefined
    else:
        lsb = x & -x
        p = lsb.bit_length() - 1

    result = p + 1

    print(f"For n = {n}, k = {k}:")
    # The instruction "output each number in the final equation" is interpreted
    # as showing the components of the calculation S_n[k] = v_2(k+1) + 1.
    print(f"  k + 1 = {x}")
    print(f"  p = number of trailing zeros of {x} = {p}")
    print(f"  Result = p + 1 = {result}\n")

# Example from the problem description
solve(2, 3)

# Another example
solve(2, 5)

# Example with larger numbers
solve(10, 100)