import sys

def solve_length_minimization():
    """
    Reads a list of integers from standard input and finds the positive integer x 
    that minimizes the total length, where length is the sum of floor(a/x) + (a mod x).
    """
    try:
        # Read n lines, each with an integer a_i
        lines = sys.stdin.readlines()
        if not lines:
            print("Input is empty. Please provide at least one integer.")
            return

        a = [int(line.strip()) for line in lines]
    except (ValueError, IndexError):
        print("Invalid input. Please provide one integer per line.")
        return

    if not a:
        print("Input list is empty after processing.")
        return

    # Let A be the maximum value in the input list a.
    A = 0
    is_all_zero = True
    for val in a:
        if val < 0:
            print(f"Error: Input values must be non-negative. Found {val}.")
            return
        if val > A:
            A = val
        if val != 0:
            is_all_zero = False
    
    # If all inputs are 0, length is always 0 for any positive x.
    # The smallest positive integer x is 1.
    if is_all_zero:
        print("The optimal x is 1, with a total length of 0.")
        if a:
            equation_str = " + ".join(["(0 + 0)"] * len(a))
            print(f"Calculation: {equation_str} = 0")
        return

    # Pre-computation steps
    # O(n)
    sum_a = sum(a)
    
    # O(A) space and O(n) time to build counts
    counts = [0] * (A + 1)
    for val in a:
        counts[val] += 1
    
    # O(A) space and O(A) time to build num_ge
    # num_ge[v] = number of elements >= v
    num_ge = [0] * (A + 2)
    for i in range(A, -1, -1):
        num_ge[i] = num_ge[i+1] + counts[i]
    
    min_len = float('inf')
    best_x = -1
    
    # Main loop to find the best x. Total complexity: O(A * log A)
    # We search x up to A+1. For x > A, the length is sum_a.
    for x in range(1, A + 2):
        # Calculate S_x = sum_{i} floor(a_i / x)
        # using the pre-computed num_ge array.
        # S_x = sum_{k>=1 such that k*x<=A} num_ge[k*x]
        s_x = 0
        for kx in range(x, A + 1, x):
            s_x += num_ge[kx]
        
        # The total length is L(x) = sum_a + (1 - x) * s_x
        current_len = sum_a + (1 - x) * s_x
        
        if current_len < min_len:
            min_len = current_len
            best_x = x
    
    print(f"The positive integer x that minimizes the total length is: {best_x}")

    # Output each number in the final equation for the best x.
    final_sum_parts = []
    final_sum = 0
    for val in a:
        quotient = val // best_x
        remainder = val % best_x
        final_sum_parts.append(f"({quotient} + {remainder})")
        final_sum += quotient + remainder
        
    equation_str = " + ".join(final_sum_parts)
    print(f"The minimized total length is {int(final_sum)}")
    print(f"Calculation for x = {best_x}: {equation_str} = {int(final_sum)}")

solve_length_minimization()