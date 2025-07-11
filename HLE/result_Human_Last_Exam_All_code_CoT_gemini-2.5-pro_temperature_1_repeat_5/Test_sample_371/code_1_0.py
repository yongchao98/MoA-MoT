def calculate_and_print_kth_element(n, k):
    """
    Calculates the k-th element of the sequence S_n using an optimal mathematical formula.

    The value of n is only to define the sequence; it's not needed for
    the optimal calculation as long as k is a valid index for S_n.
    The problem constraints imply n is large enough.

    The value of S_n[k] is given by the formula v_2(k+1) + 1, where v_2(x)
    is the number of trailing zeros in the binary representation of x.

    Args:
        n (int): The sequence number.
        k (int): The 0-based index of the element to find.
    """
    val_k_plus_1 = k + 1

    # v_2(x) can be computed efficiently using bitwise operations.
    # The expression (x & -x) isolates the lowest set bit of x.
    # The result is a power of 2, let's say 2^p.
    # The number of trailing zeros, v_2(x), is p.
    # In Python, for a number y = 2^p, its bit length is p+1.
    # So, p = y.bit_length() - 1.
    lowest_set_bit = val_k_plus_1 & -val_k_plus_1
    v2 = lowest_set_bit.bit_length() - 1

    # The final result for S_n[k] is v_2(k+1) + 1
    result = v2 + 1

    # Per the request to "output each number in the final equation",
    # we print the components of the equation: v_2(k+1) + 1 = result.
    # The numbers are v2, 1, and the result.
    print(f"{v2} + 1 = {result}")

# --- Example Usage ---
# The user can replace these values with their own.
# The problem gives an example for n=2, k=3.
# S_2 = [1, 2, 1, 3, 1, 2, 1], so S_2[3] should be 3.
n_val = 2
k_val = 3

# For the example case, the calculation is:
# k=3, k+1=4. Binary of 4 is '100'.
# v_2(4) = 2.
# Result = 2 + 1 = 3.
# The program should print "2 + 1 = 3".
calculate_and_print_kth_element(n_val, k_val)

# For a larger example, e.g., k=12345
# k+1 = 12346 = 2 * 6173. v_2(12346) = 1. Result should be 1+1=2.
# calculate_and_print_kth_element(100, 12345) # Expected output: "1 + 1 = 2"