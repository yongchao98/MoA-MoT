import math

def find_kth_element_optimal(n, k):
    """
    Calculates the k-th element of the sequence S_n using the optimal formula.
    S_n[k] = v_2(k+1) + 1, where v_2(x) is the number of trailing zeros in the
    binary representation of x.

    Args:
        n (int): The sequence number.
        k (int): The 0-based index of the element to find.

    Returns:
        int: The value of the element S_n[k].
    """
    if k < 0:
        # Assuming k is always non-negative as per sequence indexing.
        return -1 # Or raise an error
        
    val = k + 1
    
    # An efficient way to calculate v_2(val) is to find the lowest set bit
    # and take its log base 2.
    # lowest_set_bit = val & -val
    # v2 = int(math.log2(lowest_set_bit))
    # In Python 3.10+, math.cbrt() could be used as well.

    # A more direct and portable way in Python is to count trailing zeros.
    # For Python 3.10+ you can use (val).bit_count(), but that counts all set bits.
    # The `(x & -x).bit_length() - 1` trick is a great way to get v_2(x).
    if val == 0:
      return 0 # v2 of 0 is undefined, but k >= 0 so val >= 1
    
    lowest_set_bit = val & -val
    result = lowest_set_bit.bit_length() - 1 + 1

    # The final equation output format
    print(f"S_{n}[{k}] = {result}")

# Example from the problem description
n1, k1 = 2, 3
find_kth_element_optimal(n1, k1)

# Another example: S_2 = [1, 2, 1, 3, 1, 2, 1], find element at index 5
n2, k2 = 2, 5
find_kth_element_optimal(n2, k2)

# A more complex example: S_3[13]
# S_3 = [1, 2, 1, 3, 1, 2, 1, 4, 1, 2, 1, 3, 1, 2, 1]
n3, k3 = 3, 13
find_kth_element_optimal(n3, k3)
