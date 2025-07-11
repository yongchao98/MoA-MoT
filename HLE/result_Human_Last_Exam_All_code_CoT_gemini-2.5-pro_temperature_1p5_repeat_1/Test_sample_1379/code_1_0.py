def solve_tiling():
    """
    Calculates the maximal area of a rectangle that can be covered by 2x1 bricks.
    """
    N = 1000000001
    M = 1000000001

    # The total area of the rectangle is N * M.
    # Each brick has an area of 2. The total covered area must be a multiple of 2.
    # The maximal covered area can be calculated by finding the total area
    # and rounding it down to the nearest even number.
    # In integer arithmetic, this is equivalent to (N * M) // 2 * 2.
    total_area = N * M
    max_covered_area = (total_area // 2) * 2

    # The request asks to output each number in the final equation.
    print(f"Given N = {N} and M = {M}")
    print(f"The total area is {N} * {M} = {total_area}")
    print(f"The maximal coverable area with 2x1 bricks is {max_covered_area}")

solve_tiling()
