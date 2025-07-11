def find_kth_element(n, k):
    """
    Calculates the k-th element of the recursively defined sequence S_n.
    The sequence is defined as S_0 = [1] and S_n = S_{n-1} + [n+1] + S_{n-1}.

    Args:
        n (int): The sequence number, used for printing the output.
        k (int): The 0-based index of the element to find.
    """
    # The value at index k is v_2(k+1) + 1, where v_2(x) is the number of
    # trailing zeros in the binary representation of x. This can be calculated
    # efficiently for large integers in Python.

    # Let x = k + 1. The expression (x & -x) isolates the least significant
    # bit of x. For example, if x is 12 (1100), x & -x results in 4 (0100).
    # The bit_length() of a power of two, 2^p, is p+1.
    # Thus, ((k+1) & -(k+1)).bit_length() computes v_2(k+1) + 1 directly.
    
    val = k + 1
    
    # Check for invalid k, although the problem implies k >= 0.
    if val <= 0:
        result = "Invalid (k must be non-negative)"
    else:
        lsb = val & -val
        result = lsb.bit_length()

    print(f"S_{n}[{k}] = {result}")

# Main execution block with an example from the problem description
if __name__ == '__main__':
    # Example: For n=2 and k=3, calculate the 3rd element of S_2.
    # S_2 = [1, 2, 1, 3, 1, 2, 1]
    n_example = 2
    k_example = 3
    find_kth_element(n_example, k_example)