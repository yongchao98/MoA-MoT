def solve_tiling_problem():
    """
    Calculates the maximal area of an NxM rectangle that can be covered by 2x1 bricks
    and the memory usage of an equivalent C program.
    """
    # Given dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # Calculate the maximal covered area.
    # If both N and M are odd, the total area is odd, and one 1x1 square
    # cannot be covered. Otherwise, the entire area can be covered.
    if n % 2 != 0 and m % 2 != 0:
        max_covered_area = n * m - 1
    else:
        max_covered_area = n * m

    # Calculate the memory usage for an optimal C program.
    # The variables n, m, and the output would be stored as 64-bit integers (long long).
    # A 64-bit integer uses 8 bytes.
    # Memory for n: 8 bytes
    # Memory for m: 8 bytes
    # Memory for output: 8 bytes
    memory_usage_bytes = 8 + 8 + 8

    # Print the result in the format "x:o"
    # x is memory in bytes, o is the maximal covered area.
    print(f"{memory_usage_bytes}:{max_covered_area}")

solve_tiling_problem()