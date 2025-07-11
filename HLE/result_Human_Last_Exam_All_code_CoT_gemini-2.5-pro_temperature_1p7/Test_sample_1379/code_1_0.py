def solve_brick_problem():
    """
    Calculates the maximal coverable area of an N x M rectangle with 2x1 bricks
    and analyzes the memory usage of an equivalent optimal C program.
    """
    # Problem parameters
    N = 1000000001
    M = 1000000001

    # --- Part 1: Calculate the maximal area (o) ---

    # A brick's area is 2. The covered area must be even.
    # If N * M is odd (only if N and M are both odd), one 1x1 square is left.
    # For the given N and M, both are odd.
    # Python's arbitrary-precision integers handle the large multiplication easily.
    maximal_area = N * M - 1

    # --- Part 2: Calculate memory usage (x) ---

    # An optimal C program would use 64-bit integers for robustness and portability.
    # sizeof(long long) for n is 8 bytes.
    # sizeof(long long) for m is 8 bytes.
    # sizeof(unsigned long long) for the result is 8 bytes.
    memory_used_in_bytes = 8 + 8 + 8

    # --- Part 3: Output the results as specified ---

    # The prompt requires printing each number in the final equation.
    print(f"Equation for the maximal covered area:")
    print(f"{N} * {M} - 1 = {maximal_area}")
    print() # Adding a blank line for clarity

    # The final answer must be in the format x:o
    print(f"Final Answer (memory in bytes : area):")
    print(f"{memory_used_in_bytes}:{maximal_area}")

solve_brick_problem()