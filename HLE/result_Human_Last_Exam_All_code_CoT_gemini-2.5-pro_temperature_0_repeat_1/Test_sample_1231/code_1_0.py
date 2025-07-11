def solve():
    """
    This script solves the grid transformation problem by following these steps:
    1.  Define the two possible initial 8x8 grids based on the checkerboard marking pattern.
    2.  Apply the specified transformations (reflection over y=x, then 90-degree clockwise rotation)
        to each initial grid.
    3.  For each of the two resulting final grids, iterate through all 25 possible 4x4 sub-grids.
    4.  Count the number of marked cells in each 4x4 sub-grid.
    5.  Find the number of sub-grids that contain exactly 2 marked cells.
    6.  The final answer is the maximum of the counts from the two initial configurations.
    """

    def create_initial_grid(start_tl_is_marked):
        """Creates the initial 8x8 grid."""
        grid = [[0] * 8 for _ in range(8)]
        for i in range(4):
            for j in range(4):
                is_even_subgrid = (i + j) % 2 == 0
                if (start_tl_is_marked and is_even_subgrid) or \
                   (not start_tl_is_marked and not is_even_subgrid):
                    # Mark Top-Left
                    grid[2 * i][2 * j] = 1
                else:
                    # Mark Bottom-Right
                    grid[2 * i + 1][2 * j + 1] = 1
        return grid

    def transform_grid(initial_grid):
        """Applies reflection over y=x then 90-degree clockwise rotation."""
        final_grid = [[0] * 8 for _ in range(8)]
        for r in range(8):
            for c in range(8):
                # A point (col, row) reflects to (row, col).
                # Rotating (row, col) 90-deg clockwise gives (col, 7-row).
                # So, value at initial_grid[row][col] moves to final_grid[7-row][col].
                final_grid[7 - r][c] = initial_grid[r][c]
        return final_grid

    def count_subgrids_with_n_marks(grid, n):
        """Counts the number of 4x4 sub-grids with exactly n marked cells."""
        count = 0
        for r_start in range(5):
            for c_start in range(5):
                marked_cells = 0
                for r in range(r_start, r_start + 4):
                    for c in range(c_start, c_start + 4):
                        marked_cells += grid[r][c]
                if marked_cells == n:
                    count += 1
        return count

    # Case 1: Checkerboard pattern starts with a Top-Left mark.
    initial_grid_A = create_initial_grid(start_tl_is_marked=True)
    final_grid_A = transform_grid(initial_grid_A)
    count_A = count_subgrids_with_n_marks(final_grid_A, 2)

    # Case 2: Checkerboard pattern starts with a Bottom-Right mark.
    initial_grid_B = create_initial_grid(start_tl_is_marked=False)
    final_grid_B = transform_grid(initial_grid_B)
    count_B = count_subgrids_with_n_marks(final_grid_B, 2)

    result = max(count_A, count_B)

    print(f"For the first possible initial pattern, the number of 4x4 sub-grids with exactly 2 marked cells is {count_A}.")
    print(f"For the second possible initial pattern, the number of 4x4 sub-grids with exactly 2 marked cells is {count_B}.")
    print(f"The problem asks for the maximum number, which is max({count_A}, {count_B}).")
    print(f"The final answer is: {result}")

solve()