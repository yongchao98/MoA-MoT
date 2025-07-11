import sys

def get_kth_element(n, k):
    """
    Calculates the k-th element of the recursively defined sequence S_n.
    The value is determined to be v_2(k+1) + 1, where v_2(m) is the number
    of trailing zeros in the binary representation of m.
    This can be calculated efficiently using bitwise operations.
    """

    if k < 0:
        print("Error: Index k cannot be negative.")
        return

    # The maximum index for S_n is 2**(n+1) - 2.
    # We don't need to check this explicitly as the problem guarantees k is a valid index.
    
    # The value of the element is v_2(k+1) + 1.
    # Let's compute this step-by-step as requested.
    
    print(f"Given n = {n}, k = {k}")

    # Step 1: Calculate k + 1
    k_plus_1 = k + 1
    print(f"k + 1 = {k_plus_1}")
    
    # Step 2: Isolate the lowest set bit of (k+1).
    # This is equivalent to 2**v_2(k+1).
    # In Python's arbitrary precision integers, this works correctly.
    lowest_set_bit = k_plus_1 & -k_plus_1
    print(f"The lowest set bit of (k+1) is (k+1) & -(k+1) = {lowest_set_bit}")

    # Step 3: Find the exponent, which is v_2(k+1).
    # For a power of 2, say 2**p, its bit_length is p+1.
    # So, p = (2**p).bit_length() - 1
    v2 = lowest_set_bit.bit_length() - 1
    print(f"The exponent v2 = lowest_set_bit.bit_length() - 1 = {v2}")

    # Step 4: The result is v2 + 1
    result = v2 + 1
    print(f"The final result is v2 + 1 = {result}")
    
    print("-" * 20)
    print(f"The element S_{n}[{k}] is: {result}")


if __name__ == '__main__':
    # Example from the problem description
    get_kth_element(n=2, k=3)
    
    print("\n" + "="*40 + "\n")
    
    # A more complex example with larger numbers.
    # Let's find the k-th element of S_100, for k = 2^70 - 1.
    # v_2(k+1) = v_2(2^70) = 70. The result should be 70 + 1 = 71.
    # Since 71 <= 100+1, this is a valid result.
    n_large = 100
    k_large = (1 << 70) - 1  # 2**70 - 1
    get_kth_element(n_large, k_large)
