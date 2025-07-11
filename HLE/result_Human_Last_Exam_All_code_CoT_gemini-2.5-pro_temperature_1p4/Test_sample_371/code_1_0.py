def find_kth_element(n, k):
    """
    Calculates the k-th element of the sequence S_n using an optimal mathematical formula.

    The sequence is defined by S_0 = [1] and S_n = S_{n-1} + [n+1] + S_{n-1}.
    Analysis reveals that S_n[k] is equivalent to p+1, where p is the 0-indexed position
    of the least significant 1-bit in the binary representation of k+1.
    This is calculated efficiently using bitwise operations, which is the optimal approach.
    """
    print(f"Calculating the element at index k={k} for sequence S_{n}:")
    
    # Step 1: The value depends on k+1.
    val = k + 1
    
    # Step 2: The core of the calculation is a bitwise 'equation' to find the lowest set bit of k+1.
    # The operation `x & -x` isolates the lowest set bit. For example, for 12 (1100), it yields 4 (0100).
    lowest_set_bit = val & -val
    
    # Step 3: The bit length of a power of two, 2^p, is p+1. This is the final answer.
    result = lowest_set_bit.bit_length()
    
    print(f"The final equation is based on bitwise arithmetic: result = bit_length( (k+1) & -(k+1) )")
    print(f"Breaking it down:")
    print(f"k+1 = {k+1}")
    print(f"(k+1) & -(k+1) = {lowest_set_bit}")
    print(f"bit_length({lowest_set_bit}) = {result}")
    print(f"The element S_{n}[{k}] is {result}.")

# Example from the problem description: n=2, k=3. S_2[3] should be 3.
find_kth_element(n=2, k=3)

# Example where n=2, k=5. S_2 = [1, 2, 1, 3, 1, 2, 1], so S_2[5] should be 2.
find_kth_element(n=2, k=5)