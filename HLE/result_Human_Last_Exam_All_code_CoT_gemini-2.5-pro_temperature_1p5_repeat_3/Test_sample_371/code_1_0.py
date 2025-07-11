def find_kth_element(n, k):
    """
    Calculates the k-th element of the recursively defined sequence S_n.
    The optimal method is to use the formula S_n[k] = ctz(k+1) + 1,
    where ctz is the count of trailing zeros in binary representation.
    
    Args:
        n (int): The recursion depth of the sequence.
        k (int): The 0-based index of the element to find.
    """
    
    # Step 1: The element's value is determined by k, so we compute k + 1.
    x = k + 1

    # Step 2: A fast way to find the number of trailing zeros (ctz) is to isolate
    # the lowest set bit using the bitwise operation `x & -x`.
    # This works for Python's arbitrary-precision integers.
    lowest_set_bit = x & -x

    # The ctz is the base-2 logarithm of this isolated lowest set bit.
    # For a power of two, p, log2(p) can be computed as p.bit_length() - 1.
    num_trailing_zeros = lowest_set_bit.bit_length() - 1

    # Step 3: The element's value is ctz + 1.
    result = num_trailing_zeros + 1

    # Print the equation showing each number involved in the calculation, as requested.
    print(f"S_{n}[{k}] = (ctz({k} + 1)) + 1 = (ctz({x})) + 1 = {num_trailing_zeros} + 1 = {result}")

# Example from the problem description: n=2, k=3, which should result in S_2[3] = 3.
# The code will print the equation: S_2[3] = (ctz(3 + 1)) + 1 = (ctz(4)) + 1 = 2 + 1 = 3
find_kth_element(2, 3)

# Example to show another case: n=2, k=5, which should result in S_2[5] = 2.
# The code will print the equation: S_2[5] = (ctz(5 + 1)) + 1 = (ctz(6)) + 1 = 1 + 1 = 2
find_kth_element(2, 5)