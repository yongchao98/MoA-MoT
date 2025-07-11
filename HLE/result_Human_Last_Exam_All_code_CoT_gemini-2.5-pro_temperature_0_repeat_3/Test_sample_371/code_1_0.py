def solve():
    """
    This script implements the solution to find the k-th element of S_n
    and prints the result for the example case n=2, k=3.
    """

    def find_kth_element(n, k):
        """
        Calculates the k-th element of the sequence S_n.

        The value is S_n[k] = v_2(k+1) + 1, where v_2(x) is the number of
        trailing zeros in the binary representation of x. This is computed
        efficiently using bitwise operations.
        """
        # The problem uses 0-based indexing for k. The formula is based on k+1.
        val = k + 1

        # This case should not happen based on problem constraints (k>=0).
        if val <= 0:
            raise ValueError("k must be non-negative.")
            
        # The bitwise operation (val & -val) isolates the least significant bit (LSB).
        # For a number 'p' that is a power of two, p.bit_length() is log2(p) + 1.
        # The number of trailing zeros is v_2(val) = log2(LSB) = LSB.bit_length() - 1.
        # The final result is v_2(val) + 1, which simplifies to LSB.bit_length().
        lsb = val & -val
        result = lsb.bit_length()
        return result

    # Example from the prompt
    n_example = 2
    k_example = 3

    # Calculate the result for the example
    result_example = find_kth_element(n_example, k_example)

    # Print the final equation as requested
    print(f"S_{n_example}[{k_example}] = {result_example}")

solve()