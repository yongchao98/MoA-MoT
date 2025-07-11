import math

def solve():
    """
    Finds the positive integer x that minimizes the total length for a list of integers.
    The input is read from stdin. First line is n, followed by n lines each with an integer.
    """
    import sys
    
    try:
        # Fast I/O
        input = sys.stdin.readline
        lines = sys.stdin.readlines()
        
        # It's common to have n on the first line, then n numbers.
        # This code assumes all lines are numbers `a_i`.
        # If the first line is n, uncomment the line below.
        # n = int(lines[0])
        # a = [int(line) for line in lines[1:]]
        a = [int(line) for line in lines if line.strip()]
        if not a:
            print("Input is empty.")
            return

    except (IOError, ValueError) as e:
        print(f"Error reading input: {e}")
        return

    # Find the maximum value in a, which we denote as A in complexity analysis.
    max_a = 0
    for val in a:
        if val > max_a:
            max_a = val

    candidates = {1}
    # Generate candidate values for x
    # For each a_i, the interesting values for x are around divisors of a_i.
    # The value of floor(a_i/x) changes only O(sqrt(a_i)) times.
    # The critical points are of the form floor(a_i/k).
    # We generate these by iterating k from 1 to sqrt(a_i).
    for val in a:
        limit = int(math.sqrt(val))
        for k in range(1, limit + 1):
            candidates.add(k)
            # Add the quotient floor(a_i / k) as a candidate
            candidates.add(val // k)

    min_total_length = float('inf')
    best_x = -1

    # Sort candidates to check them in an orderly fashion, though not strictly necessary.
    sorted_candidates = sorted(list(candidates))

    for x in sorted_candidates:
        if x == 0:  # x must be a positive integer
            continue
            
        current_total_length = 0
        is_valid_candidate = True
        for val in a:
            # The calculation is q + r, where q = floor(a/x) and r = a mod x.
            current_total_length += val // x + val % x
            
        if current_total_length < min_total_length:
            min_total_length = current_total_length
            best_x = x

    # Final result construction is not required by the prompt,
    # but showing the process for finding the min sum.
    # We construct the final equation with the best_x found.
    # For demonstration, we'll print the result in a readable format.
    
    final_sum_str_parts = []
    final_sum_calc = 0
    for val in a:
        q = val // best_x
        r = val % best_x
        final_sum_str_parts.append(f"({q} + {r})")
        final_sum_calc += q + r
    
    # print(f"The integer x that minimizes the total length is: {best_x}")
    # print(f"The minimized total length is: {min_total_length}")
    # print("Calculation:")
    # print(f"sum = {' + '.join(final_sum_str_parts)} = {final_sum_calc}")

solve()
# The problem asks for the time complexity of the best algorithm.
# Based on our analysis, a well-implemented algorithm has a complexity
# related to the number of candidates and the size of the input list.
# The most practical bound considering common problem constraints is O(n*sqrt(A)).

<<<O(n*sqrt(A))>>>