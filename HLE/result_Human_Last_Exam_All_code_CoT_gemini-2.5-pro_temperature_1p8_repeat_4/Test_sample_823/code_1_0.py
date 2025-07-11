import math

def calculate_induced_matching_size_in_grid(k):
    """
    Calculates the size of an induced matching in a k-by-k grid.

    The construction places horizontal edges starting at (3i, 3j) for i,j = 0, 1, 2, ...
    This spacing ensures no two edges in the matching are adjacent.
    """
    if k <= 0:
        print("Grid size k must be a positive integer.")
        return

    # We space the matching edges by 3 units (horizontally and vertically)
    # to ensure the matching is induced.
    # The number of edges we can place in one dimension is floor(k/3).
    m = k // 3

    # The total size of the induced matching is m * m.
    size = m * m

    print(f"For a {k}x{k} grid, we can construct an induced matching.")
    print(f"The number of matching edges along one dimension is m = {k} // 3 = {m}.")
    print(f"The total size of the induced matching is {m} * {m} = {size}.")

# Example for a 20x20 grid
calculate_induced_matching_size_in_grid(20)
