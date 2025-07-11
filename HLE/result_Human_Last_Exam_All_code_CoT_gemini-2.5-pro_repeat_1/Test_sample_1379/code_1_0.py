def solve_brick_problem():
    """
    Calculates the maximal area of a rectangle that can be covered by 2x1 bricks
    and the memory required for the variables in a corresponding C program.
    """
    # Given dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # In a C program, to store numbers of this magnitude (n, m, and their product),
    # we would need to use a 64-bit integer type, such as 'long long'.
    # A 'long long' variable typically uses 8 bytes of memory.
    # The problem asks for the memory used by variables for n, m, and the output.
    # So, memory_x = sizeof(n) + sizeof(m) + sizeof(output).
    size_of_long_long_in_bytes = 8
    memory_x = 3 * size_of_long_long_in_bytes

    # The area of a single brick is 2. The total covered area must be a multiple of 2.
    # The maximal area is the total area of the rectangle, rounded down to the nearest
    # even number. This can be calculated using integer division.
    # Python's integers can handle arbitrarily large numbers, so overflow is not an issue.
    total_area = n * m
    maximal_area_o = (total_area // 2) * 2

    # Print the result in the specified format "x:o"
    print(f"{memory_x}:{maximal_area_o}")

solve_brick_problem()