import sys

def solve():
    """
    Calculates the k-th element of the sequence S_n.
    
    The function takes n and k as inputs from command-line arguments,
    or uses large default values if no arguments are provided.
    """
    
    # Use command-line arguments if provided, otherwise use large default values.
    # The problem statement implies n and k are extremely large.
    # n is conceptually infinite for this calculation, but we'll use a placeholder.
    if len(sys.argv) == 3:
        try:
            n = int(sys.argv[1])
            k = int(sys.argv[2])
        except ValueError:
            print("Please provide valid integer values for n and k.")
            return
    else:
        # Using default large values as per the problem description "n, k > INT_MAX/2"
        # INT_MAX is 2**31 - 1, so INT_MAX/2 is roughly 10**9.
        # Let's use a k value around 10^18 (which is > 2^31) for demonstration.
        n = 10**18 
        k = 10**18 + 3 # For k+1 = 10^18 + 4

    # The value of S_n[k] is independent of n (assuming n is large enough
    # for k to be a valid index) and follows the formula: v_2(k+1) + 1.
    # v_2(x) is the number of trailing zeros in the binary representation of x.

    if k < 0:
        print("Error: Index k cannot be negative.")
        return

    # Let x = k + 1
    x = k + 1

    # To find p = v_2(x), we can use a bitwise trick.
    # 'x & -x' isolates the lowest set bit of x. This results in a power of 2, 2^p.
    # For example, if x = 12 (binary '1100'), lowest_set_bit is 4 (binary '100').
    if x == 0:
        # This case occurs if k = -1, which is invalid.
        # The number of trailing zeros is undefined. We handle it for robustness.
        print("Error: Cannot calculate for k = -1 as k+1 is zero.")
        return
        
    lowest_set_bit = x & -x

    # Now, we find the exponent p from the value 2^p.
    # For a power of two y = 2^p, its bit length is p+1.
    # So, p = y.bit_length() - 1.
    p = lowest_set_bit.bit_length() - 1

    # The final result is p + 1.
    result = p + 1
    
    # Per instructions, "output each number in the final equation"
    # The final equation is S_n[k] = p + 1 = result
    print(f"S_{n}[{k}] = {p} + 1 = {result}")

solve()