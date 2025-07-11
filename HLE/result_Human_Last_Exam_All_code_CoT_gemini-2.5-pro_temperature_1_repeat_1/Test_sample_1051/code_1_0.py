import sys

def solve():
    """
    Finds the positive integer x that minimizes the total length.
    """
    # Read all lines from standard input
    try:
        lines = sys.stdin.readlines()
        if not lines:
            print("Input is empty.")
            return
        a_values = [int(line.strip()) for line in lines]
    except (ValueError, IndexError):
        print("Invalid input: all lines must contain a single integer.")
        return

    if not a_values:
        print("No numbers provided.")
        return

    A = 0
    for val in a_values:
        if val > A:
            A = val
    
    if A == 0:
        # If all a_i are 0, any x works, length is 0. We can pick x=1.
        # The equation will be a sum of zeros.
        equation_parts = ["0"] * len(a_values)
        print(" + ".join(equation_parts) + " = 0")
        return

    counts = [0] * (A + 2)
    total_sum_S = 0
    for val in a_values:
        counts[val] += 1
        total_sum_S += val

    g = [0] * (A + 2)
    for v in range(A, 0, -1):
        g[v] = g[v + 1] + counts[v]

    F = [0] * (A + 2)
    for x in range(1, A + 1):
        for kx in range(x, A + 1, x):
            F[x] += g[kx]

    min_total_length = -1
    best_x = -1

    # Check x from 1 to A
    for x in range(1, A + 1):
        total_length = total_sum_S + (1 - x) * F[x]
        if best_x == -1 or total_length < min_total_length:
            min_total_length = total_length
            best_x = x

    # Check x > A, for example x = A + 1. For such x, F(x)=0.
    total_length_A_plus_1 = total_sum_S
    if best_x == -1 or total_length_A_plus_1 < min_total_length:
        min_total_length = total_length_A_plus_1
        best_x = A + 1
        
    # Output the final equation
    equation_parts = []
    final_sum = 0
    for val in a_values:
        length_i = val // best_x + val % best_x
        equation_parts.append(str(length_i))
        final_sum += length_i
    
    print(" + ".join(equation_parts) + " = " + str(final_sum))

solve()
