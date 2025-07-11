def find_kth_element(n, k):
    """
    Calculates the k-th element of the sequence S_n using an optimal method.

    The recursive definition of the sequence is:
    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1}

    Analysis reveals that the value of the k-th element, S_n[k], follows a
    pattern related to the binary representation of the index. Specifically,
    S_n[k] is equal to the 1-based position of the least significant bit (LSB)
    of the number (k+1).

    This value can be calculated as v_2(k+1) + 1, where v_2(x) is the 2-adic
    valuation of x (i.e., the number of trailing zeros in the binary form of x).

    The implementation uses a fast bitwise trick to find this value.
    The expression (x & -x) isolates the lowest set bit of x. For a number x,
    this results in 2**v_2(x). The bit_length() of 2**p is p+1.
    Therefore, ((k+1) & -(k+1)).bit_length() computes v_2(k+1) + 1 directly.

    Args:
        n (int): The sequence number. Per the problem, this can be a very large integer.
                 It is not used in the calculation but defines the valid range for k.
        k (int): The 0-based index of the element to find. Per the problem, this
                 can be a very large integer.

    Returns:
        int: The value of the k-th element in S_n.
    """

    # The value of n is only needed to validate the index k, i.e., k < 2**(n+1) - 1.
    # The formula itself does not depend on n.

    # Let x be k+1
    x = k + 1

    # Calculate the result using the optimal bitwise method.
    result = (x & -x).bit_length()

    return result

# Example from the problem description: n = 2, k = 3.
# The sequence S_2 is [1, 2, 1, 3, 1, 2, 1].
# The element at index 3 is 3.
n_example = 2
k_example = 3
value = find_kth_element(n_example, k_example)

# The prompt requires printing the result.
# The instruction "output each number in the final equation!" is interpreted
# as printing the final computed value.
print(value)
