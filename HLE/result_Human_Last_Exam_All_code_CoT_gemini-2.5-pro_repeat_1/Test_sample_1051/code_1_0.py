import sys

def solve():
    """
    Finds the optimal x that minimizes the total length and prints the result.
    """
    # Read all lines from standard input
    lines = sys.stdin.readlines()
    if not lines:
        # Handle empty input
        print("Input is empty.")
        return
        
    try:
        a = [int(line.strip()) for line in lines]
    except ValueError:
        print("Invalid input: all lines must contain integers.")
        return

    if not a:
        print("Input list is empty after processing.")
        return

    A = 0
    for val in a:
        if val > A:
            A = val
    
    if A == 0:
        # All a_i are 0, any x>=1 gives length 0. We choose 1.
        best_x = 1
    else:
        # Precomputation for the N(y) function
        # N[y] = count of a_i >= y
        counts = [0] * (A + 1)
        for val in a:
            if val > 0:
                counts[val] += 1
        
        N = [0] * (A + 2)
        for y in range(A, 0, -1):
            N[y] = N[y+1] + counts[y]

        best_x = 1
        max_g = 0.0

        # Iterate through all possible x from 1 to A
        # This is the O(A log A) part
        for x in range(1, A + 1):
            if x == 1:
                # g(1) = 0
                continue
            
            # Calculate F(x) = sum(floor(a_i/x))
            # F(x) = sum_{k=1..A/x} N(k*x)
            F_x = 0
            # k*x <= A
            for k in range(1, A // x + 1):
                F_x += N[k * x]
            
            g_x = (x - 1) * F_x
            
            if g_x > max_g:
                max_g = g_x
                best_x = x

    # Output the result in the requested format
    total_length = 0
    equation_parts = []
    for val in a:
        q = val // best_x
        r = val % best_x
        length_i = q + r
        total_length += length_i
        equation_parts.append(f"({q} + {r})")
        
    equation_str = " + ".join(equation_parts)
    print(f"{equation_str} = {total_length}")

solve()
