def find_kth_element(n, k):
    """
    Calculates the k-th element of the sequence S_n using an optimal method.
    The value S_n[k] is determined by the number of trailing zeros in the
    binary representation of k+1.

    Args:
        n (int): The index of the sequence (note: the result is independent of n).
        k (int): The index of the element to find (0-based).

    Returns:
        The value of the k-th element in S_n.
    """
    # The value we need to analyze is k+1.
    x = k + 1

    # For any positive integer x, the expression (x & -x) in two's complement
    # arithmetic isolates the lowest set bit. For example, if x = 12 (..001100),
    # -x is (..110100), and (x & -x) is (..000100), which is 4.
    # The result is 2^p, where p is the number of trailing zeros.
    lowest_set_bit = x & -x

    # To find p from 2^p, we can use the bit_length() method.
    # (2**p).bit_length() returns p+1.
    num_trailing_zeros = lowest_set_bit.bit_length() - 1

    # The value of the element in the sequence is the number of trailing zeros plus one.
    result = num_trailing_zeros + 1

    return result

# Example from the problem description: n = 2, k = 3.
# S_2 = [1, 2, 1, 3, 1, 2, 1]. The element at index 3 is 3.
n_example = 2
k_example = 3
result_example = find_kth_element(n_example, k_example)

# As requested, printing the final equation with each number.
print(f"S_{n_example}[{k_example}] = {result_example}")