import sys

def get_kth_element(n, k):
    """
    Calculates the k-th element of the recursively defined sequence S_n.

    S_0 = [1]
    S_n = S_{n-1} + [n+1] + S_{n-1} for n > 0

    The optimal solution relies on a pattern found in the sequence's structure.
    The value S_n[k] is independent of n (for sufficiently large n) and is
    equal to v+1, where v is the number of trailing zeros in the binary
    representation of the 1-based index, k+1.
    """
    if k < 0:
        print("Error: Index k cannot be negative.")
        return

    # The value we are looking for is `v+1`, where `v` is the number of
    # trailing zeros in the binary representation of `k+1`.
    
    # Let's take the example n=2, k=3: S_2 = [1, 2, 1, 3, 1, 2, 1] -> S_2[3] = 3.
    # 1. k = 3, so k+1 = 4.
    # 2. Binary of 4 is '100'. It has v=2 trailing zeros.
    # 3. The value is v+1 = 2+1 = 3. This is correct.
    
    # We can compute this value efficiently using bitwise operations.
    
    # `k+1` converts the 0-based index to a 1-based index.
    target = k + 1
    
    # The expression `target & -target` isolates the lowest set bit.
    # For a number equal to 2^v * (something odd), this expression yields 2^v.
    lowest_set_bit = target & -target
    
    # For a power of two, 2^v, its bit_length() is v+1.
    # This value, v+1, is exactly the result we need.
    result = lowest_set_bit.bit_length()

    # The prompt asks to output each number in the final equation.
    # The final relation can be expressed as S_n[k] = result.
    # Thus, we print the values for n, k, and the calculated result.
    print(f"S_{n}[{k}] = {result}")


if __name__ == '__main__':
    # Default example from the problem description: n=2, k=3.
    n_val = 2
    k_val = 3

    # The problem can be run with command-line arguments.
    # For example: python your_script.py 100 50
    if len(sys.argv) == 3:
        try:
            n_val = int(sys.argv[1])
            k_val = int(sys.argv[2])
        except ValueError:
            print("Invalid arguments. Please provide two integers for n and k.")
            print(f"Running with default values n={n_val}, k={k_val}.")
            
    get_kth_element(n_val, k_val)