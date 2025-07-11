def solve():
    """
    Calculates the maximum possible number of digits in N based on the problem's rules.
    """
    
    # L(k) is the maximum length for k distinct digits.
    # We start with k=1. The longest sequence is just a single digit, e.g., "1".
    # The sequence "11" is invalid because the digit 1 appears twice.
    k = 1
    length = 1
    print(f"For k = {k} distinct digit:")
    print(f"The maximum length is {length}.")
    print("-" * 30)

    # Now, we calculate the lengths for k = 2, 3, 4, 5 using the recurrence L(k) = 2*L(k-1) + 1
    for k in range(2, 6):
        prev_length = length
        length = 2 * prev_length + 1
        print(f"For k = {k} distinct digits:")
        print(f"The maximum length follows the equation: 2 * {prev_length} + 1 = {length}")
        print("-" * 30)
    
    print(f"The maximum possible number of digits in N is the length for k=5, which is {length}.")

solve()