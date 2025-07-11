def solve_brick_problem():
    """
    Calculates the maximal area of a rectangle that can be covered by 2x1 bricks
    and the memory required for a C program to store the variables.
    """

    # The dimensions of the rectangle.
    # Python's integers handle arbitrary size, so no overflow will occur.
    n = 1000000001
    m = 1000000001

    # The area of one brick is 2. The total area covered must be even.
    # The total area of the rectangle is n * m.
    # The maximum coverable area is floor(total_area / 2) * 2.
    # In integer arithmetic, this is equivalent to (n * m) // 2 * 2.
    max_area = (n * m) // 2 * 2

    # The problem asks for the memory used by an optimal C program.
    # On a typical 64-bit system, the values for n, m, and the result
    # fit within the 'long long' data type, which uses 8 bytes.
    # Memory for n: 8 bytes
    # Memory for m: 8 bytes
    # Memory for the output: 8 bytes
    # Total memory: 8 + 8 + 8 = 24 bytes.
    memory_used_in_bytes = 24

    # The final answer should be in the format x:o,
    # where x is memory in bytes and o is the maximal area.
    print(f"{memory_used_in_bytes}:{max_area}")

solve_brick_problem()