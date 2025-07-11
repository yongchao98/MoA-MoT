def solve_tiling_problem():
    """
    Calculates the maximal covered area of an N x M rectangle with 2x1 bricks
    and the memory usage of an equivalent C program.
    """

    # The dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # The product N * M will be odd since both N and M are odd.
    # The maximal area that can be covered by 2x1 bricks is N * M - 1.
    # In integer arithmetic, (n * m // 2) * 2 achieves this correctly.
    maximal_area = (n * m // 2) * 2

    # In a C program, N, M, and the result (N*M) would require
    # a 64-bit integer type (long long) to avoid overflow.
    # The size of a long long is 8 bytes on most systems.
    # We need three such variables: one for n, one for m, and one for the output.
    size_of_long_long = 8  # in bytes
    num_variables = 3
    memory_used = size_of_long_long * num_variables

    # The problem requires the output in the format x:o
    # x = memory used in bytes
    # o = the calculated maximal area
    print(f"{memory_used}:{maximal_area}")

solve_tiling_problem()