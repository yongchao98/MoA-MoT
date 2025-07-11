def solve_kth_element(n, k):
    """
    Calculates the k-th element of the sequence S_n, where S_0 = [1]
    and S_n = S_{n-1} + [n+1] + S_{n-1}.
    The k-th element (0-indexed) of S_n follows the pattern: S_n[k] = tz(k+1) + 1,
    where tz(x) is the number of trailing zeros in the binary representation of x.

    Args:
        n: The sequence number (a non-negative integer).
        k: The index of the element to find (a non-negative integer).
    """

    # The value of n is primarily for context, ensuring k is a valid index.
    # The actual calculation depends only on k.
    
    # We need to find the number of trailing zeros of k+1.
    val = k + 1
    
    # A bitwise trick `x & -x` can isolate the least significant bit (LSB) of x.
    # For val > 0, lsb will be a power of two, 2**p, where p is the number of trailing zeros.
    if val <= 0:
        # According to the problem, k is a 0-based index, so k>=0 and val>=1.
        # This case handles unexpected inputs.
        print("Error: k must be non-negative.")
        return

    lsb = val & -val
    
    # The bit_length of 2**p is p+1.
    # So, the number of trailing zeros, p, is lsb.bit_length() - 1.
    trailing_zeros = lsb.bit_length() - 1
    
    # The value of the element is tz(k+1) + 1.
    result = trailing_zeros + 1
    
    # Print the final equation with all numbers.
    print(f"S_{n}[{k}] = tz({k}+1) + 1 = tz({val}) + 1 = {trailing_zeros} + 1 = {result}")

# Example from the problem description: n = 2, k = 3
n_example = 2
k_example = 3
solve_kth_element(n_example, k_example)