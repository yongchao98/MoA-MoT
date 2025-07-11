import math
from collections import Counter

def solve():
    """
    Finds the positive integer x that minimizes the total length for a list of integers.
    The input is read from standard input, with the first line being the number of integers n,
    followed by n lines each containing an integer a_i.
    """
    try:
        # Read input from the user
        # In an interactive environment, you would type the numbers and then press Ctrl+D (on Linux/macOS) or Ctrl+Z+Enter (on Windows)
        import sys
        input_lines = sys.stdin.readlines()
        if not input_lines:
            # Default example if no input is provided
            print("No input detected, running with a default example: a = [13, 23, 17, 9, 29]")
            a = [13, 23, 17, 9, 29]
        else:
            # First line is often n, but we can just parse all lines as integers
            a = [int(line.strip()) for line in input_lines if line.strip().isdigit()]
            if not a and len(input_lines) > 0 and input_lines[0].strip().isdigit():
                 n = int(input_lines[0].strip())
                 # assume a_i are on the next n lines, already read. Let's re-parse
                 a = [int(line.strip()) for line in input_lines[1:] if line.strip().isdigit()]


    except (ValueError, IndexError):
        # Fallback for incorrect input format or execution environment
        print("Invalid input, running with a default example: a = [13, 23, 17, 9, 29]")
        a = [13, 23, 17, 9, 29]

    if not a:
        print("The list of integers is empty.")
        return

    counts = Counter(a)
    unique_a = list(counts.keys())
    
    # --- Candidate Generation ---
    candidates = {1}
    for v in unique_a:
        if v == 0:
            continue
        sqrt_v = int(math.sqrt(v))
        for k in range(1, sqrt_v + 1):
            candidates.add(k)
            candidates.add(v // k)
            if v // k > 1:
                candidates.add((v // k) - 1)
    
    candidates.discard(0)

    # --- Find Best x ---
    min_total_length = float('inf')
    best_x = -1
    
    # Pre-calculate Sum(a_i)
    total_sum_a = sum(val * count for val, count in counts.items())

    for x in sorted(list(candidates)):
        q_sum = 0
        for v, count in counts.items():
            q_sum += count * (v // x)
        
        current_total_length = total_sum_a + q_sum * (1 - x)
        
        if current_total_length < min_total_length:
            min_total_length = current_total_length
            best_x = x

    # --- Output Results ---
    print(f"The optimal value for x is: {best_x}")
    print(f"The minimized total length is: {min_total_length}")

    # For display purposes, calculate the breakdown of the sum
    final_sum_parts = []
    # Using the original list 'a' to preserve order for the equation
    for val in sorted(a):
        q = val // best_x
        r = val % best_x
        final_sum_parts.append(f"({q} + {r})")

    equation_str = " + ".join(final_sum_parts)
    print(f"\nThe equation with x = {best_x} is:")
    print(f"Total Length = {equation_str} = {min_total_length}")


solve()