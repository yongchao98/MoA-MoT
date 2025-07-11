def solve_brick_problem():
    """
    Calculates memory usage for a C program and the max tileable area.
    """
    
    # 1. Define the dimensions of the rectangle
    N = 1000000001
    M = 1000000001

    # 2. Calculate the memory 'x' that an optimal C program would use.
    # The variables are n, m, and the output.
    # N, M, and their product (N*M) require a 64-bit integer type (e.g., unsigned long long in C).
    # A 64-bit integer typically occupies 8 bytes.
    # memory for n (8 bytes) + memory for m (8 bytes) + memory for output (8 bytes)
    memory_in_bytes = 8 + 8 + 8

    # 3. Calculate the maximal coverable area 'o'.
    # The formula using integer arithmetic is (N * M) / 2 * 2.
    # This finds the largest even number less than or equal to N * M.
    maximal_area = (N * M) // 2 * 2

    # 4. Print the final result in the format "x:o".
    # The equation solved is (1000000001 * 1000000001) // 2 * 2
    print(f"{memory_in_bytes}:{maximal_area}")

solve_brick_problem()