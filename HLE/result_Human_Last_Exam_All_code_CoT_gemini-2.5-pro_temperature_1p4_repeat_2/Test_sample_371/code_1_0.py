def find_kth_element(n, k):
    """
    Calculates the k-th element of the recursively defined sequence S_n.
    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1} for n > 0

    Args:
        n: The sequence number (an integer, can be very large).
        k: The 0-based index of the element to find (an integer, can be very large).

    Returns:
        The integer value of the k-th element in S_n.
    """

    # The value S_n[k] is independent of n (for valid k) and follows
    # the pattern of the ruler function.
    # It is equal to the 1-based index of the least significant '1' bit
    # in the binary representation of (k+1).

    # Let x = k + 1.
    x = k + 1

    # The value of the least significant bit can be isolated using the
    # bitwise operation `x & -x`. For example, if x is 12 (0b1100),
    # `x & -x` results in 4 (0b100).
    lsb_power_of_2 = x & -x

    # The result we need is the position of this bit. For a number that is
    # a power of 2 (like lsb_power_of_2), its 1-based bit position is
    # equal to its bit_length(). For example, 4 (0b100) is the 3rd bit,
    # and 4.bit_length() returns 3.
    return lsb_power_of_2.bit_length()

# Example from the problem description: n=2, k=3
# We expect S_2[3] = 3.
# Let's verify with the function:
# k=3 -> k+1 = 4.
# 4 in binary is 100.
# lsb_power_of_2 = 4 & -4 = 4.
# 4.bit_length() is 3. The function returns 3.
n_example = 2
k_example = 3
result = find_kth_element(n_example, k_example)
print(result)
