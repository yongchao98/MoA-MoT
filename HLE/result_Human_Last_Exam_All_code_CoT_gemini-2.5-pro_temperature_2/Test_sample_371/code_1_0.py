import sys

def get_kth_element():
    """
    Reads n and k from standard input and calculates the k-th element of S_n.

    The sequence S_n is defined as:
    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1} for n > 0

    The value of the k-th element S_n[k] is found to be ctz(k+1) + 1,
    where ctz is the count of trailing zeros in the binary representation of a number.
    
    This function expects two space-separated integers (n and k) from stdin.
    """
    try:
        # Read n and k from input. They can be very large.
        n_str, k_str = sys.stdin.readline().split()
        n = int(n_str)
        k = int(k_str)
    except (ValueError, IndexError):
        # Provide default values if input is missing or invalid, as per example
        n = 2
        k = 3
        print(f"Using default values n={n}, k={k} since input was not provided.", file=sys.stderr)


    # The value of S_n[k] depends only on k. The value of n defines the sequence
    # and constrains the maximum value of k, but isn't needed for the calculation
    # itself, assuming k is a valid index.
    
    # We need to calculate ctz(k + 1). We add 1 to k for the calculation.
    val = k + 1
    
    # Now, we count the number of trailing zeros (ctz).
    # A simple way is to count how many times we can divide by 2 (right shift)
    # before the number becomes odd.
    count = 0
    
    # Since k is a non-negative index, k+1 is always >= 1.
    # The loop condition `(val & 1) == 0` checks if the number is even.
    while (val & 1) == 0:
        val >>= 1
        count += 1
        
    # The value of the element is ctz(k+1) + 1.
    result = count + 1
    
    print(result)

# Execute the function
get_kth_element()
