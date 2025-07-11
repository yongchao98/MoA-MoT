def solve_kth_element():
    """
    This script defines and executes the optimal function to find S_n[k].
    """

    def find_kth_element_optimal(n, k):
        """
        Calculates the k-th element of the sequence S_n using the optimal method.
        
        The value S_n[k] is equivalent to the ruler function, which is the
        number of trailing zeros in the binary representation of the 1-based index (k+1),
        plus one.
        
        Args:
            n: The sequence number (assumed to be very large).
            k: The 0-based index of the element to find (assumed to be very large).
            
        Returns:
            The integer value of the k-th element.
        """
        # The problem states k is a 0-based index.
        # The pattern relates to the 1-based index, j.
        j = k + 1
        
        # The number of trailing zeros (p) in j can be found efficiently.
        # The expression (j & -j) isolates the lowest set bit of j.
        # For j > 0, the value of this expression is 2**p.
        # The bit_length() of 2**p is p + 1.
        # Therefore, p = (j & -j).bit_length() - 1.
        if j <= 0:
            # This case should not happen given k >= 0.
            p = -1 
        else:
            p = (j & -j).bit_length() - 1
            
        result = p + 1
        
        # As requested, output the components of the final calculation.
        # We use string representations for potentially huge n and k.
        n_str = str(n)
        k_str = str(k)
        if len(n_str) > 10: n_str = "n"
        if len(k_str) > 10: k_str = "k"
            
        print(f"S_{n_str}[{k_str}] = ctz({k} + 1) + 1 = ctz({j}) + 1 = {p} + 1 = {result}")
        
        return result

    # Example from the problem description: n = 2, k = 3 -> S_2[3] = 3
    print("Example calculation for n=2, k=3:")
    find_kth_element_optimal(2, 3)

# Execute the solution
solve_kth_element()