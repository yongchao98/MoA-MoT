def get_kth_element_optimal(n, k):
    """
    Calculates the k-th element of the sequence S_n using an optimal bitwise method.

    The value of the k-th element (0-indexed) corresponds to the 1-indexed position
    of the least significant bit (LSB) in the binary representation of (k+1).

    Args:
        n (int): The level of the sequence. Used for the output string.
        k (int): The 0-based index of the element to find.

    Returns:
        int: The value of the k-th element.
    """
    if k < 0:
        raise ValueError("Index k cannot be negative.")

    # Python's integers handle arbitrary size, so k can be very large.
    k_plus_1 = k + 1

    # Isolate the lowest set bit. For an integer x, (x & -x) returns
    # a power of 2 corresponding to the LSB of x.
    # For example, if k+1 = 6 (binary 110), lowest_set_bit is 2 (binary 10).
    lowest_set_bit = k_plus_1 & -k_plus_1

    # The position of this bit (0-indexed) can be found by taking log base 2.
    # In Python, for a power of two `v = 2**p`, `v.bit_length() - 1` efficiently computes `p`.
    # This value `p` is the 0-indexed position of the LSB.
    p = lowest_set_bit.bit_length() - 1

    # The problem's sequence value is p + 1.
    result = p + 1
    
    # Per the instructions, output the numbers in the final equation.
    # We will format this as S_n[k] = result.
    print(f"S_{n}[{k}] = {result}")
    
    return result

# Example from the problem description to demonstrate the function
# n = 2, k = 3, S_2[3] should be 3.
# k+1 = 4 (binary 100). LSB is at position 2 (0-indexed). Result = 2+1=3.
get_kth_element_optimal(2, 3)

# Another example: n=2, k=5. S_2 = [1, 2, 1, 3, 1, 2, 1]. S_2[5] should be 2.
# k+1 = 6 (binary 110). LSB is at position 1 (0-indexed). Result = 1+1=2.
get_kth_element_optimal(2, 5)

# Example with a larger n and k, assuming they are within standard integer types for demonstration.
# In the actual problem, these would be huge numbers.
get_kth_element_optimal(10, 100)