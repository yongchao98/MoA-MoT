import math

def solve_for_kth_element():
    """
    This function calculates the k-th element of a recursively defined sequence S_n.
    The sequence is defined by S_0 = [1] and S_n = S_{n-1} + [n+1] + S_{n-1}.

    The optimal solution finds that the value of S_n[k] is independent of n
    and is equal to the number of trailing zeros in the binary representation
    of (k+1), plus one. This is denoted as ctz(k+1) + 1.
    """
    # Example values from the problem description
    n = 2
    k = 3

    # The equation to find the result is: S_n[k] = ctz(k + 1) + 1.
    # We will output the numbers involved in this calculation.

    # 1. The first number is k+1
    k_plus_1 = k + 1

    # 2. The second number is ctz(k+1). This can be calculated efficiently
    #    using bitwise operations. `x & -x` isolates the lowest set bit of x,
    #    which is a power of 2 (2^ctz). The bit_length of a power of two minus 1
    #    gives its exponent.
    if k_plus_1 <= 0:
        # This case is not expected for k >= 0.
        ctz = -1 
    else:
        lowest_set_bit = k_plus_1 & -k_plus_1
        ctz = lowest_set_bit.bit_length() - 1
    
    # 3. The final number is the result.
    result = ctz + 1

    # As requested, output each number in the final equation.
    # The final equation is: result = ctz + 1
    # The numbers involved are n, k, k+1, ctz, and the final result.
    print(f"n = {n}")
    print(f"k = {k}")
    print(f"k + 1 = {k_plus_1}")
    print(f"ctz(k + 1) = {ctz}")
    print(f"Final result = {result}")

solve_for_kth_element()