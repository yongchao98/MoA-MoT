import math

def solve():
    """
    Calculates the k-th element of the sequence S_n based on the optimal
    method and prints the result as an equation.
    
    The sequence S_n has a property where its k-th element (0-indexed) is
    equal to the number of trailing zeros in the binary representation of k+1, plus 1.
    S_n[k] = ctz(k+1) + 1
    
    This avoids generating the sequence, which would be infeasible for large n.
    """
    # Example values as given in the problem description.
    n = 2
    k = 3

    # The value is independent of n, as long as k is a valid index for S_n.
    x = k + 1

    # An efficient way to count trailing zeros (ctz) for a positive integer x.
    # Since Python 3.10, math.ctz(x) can be used. For compatibility, we use
    # a bit manipulation trick.
    # x & -x isolates the least significant bit. The number of trailing zeros
    # is the base-2 logarithm of this value.
    if x <= 0:
        # This case applies if k is -1 or less. Assuming k is a valid non-negative index.
        num_trailing_zeros = -1 # Or handle as an error
    else:
        lowest_set_bit = x & -x
        num_trailing_zeros = lowest_set_bit.bit_length() - 1

    result = num_trailing_zeros + 1
    
    # The prompt requests to "output each number in the final equation".
    # The final calculation is `result = num_trailing_zeros + 1`.
    # We print the numbers involved in this equation.
    print(f"{num_trailing_zeros} + 1 = {result}")

solve()