import numpy as np

def solve_grid_problem():
    """
    This function solves the grid transformation problem by:
    1. Creating an initial 8x8 grid with marked cells based on the specified pattern.
    2. Applying reflection and rotation transformations to the marked cells.
    3. Counting how many 4x4 sub-grids in the final configuration have exactly two marked cells.
    """
    N = 8
    subgrid_size = 4

    # Step 1: Generate the initial marked cells
    # The 8x8 grid consists of 4x4 = 16 sub-grids of size 2x2.
    # We can index these sub-grids by (I, J) where I, J are from 0 to 3.
    # The checkerboard pattern places marks at (i,j) or (i+1, j+1)
    # depending on whether I+J is even or odd.
    initial_marks = []
    for I in range(N // 2):
        for J in range(N // 2):
            # Top-left corner of the 2x2 sub-grid
            i, j = 2 * I, 2 * J
            
            if (I + J) % 2 == 0:
                # Mark the top-left corner of the 2x2 sub-grid
                initial_marks.append((i, j))
            else:
                # Mark the bottom-right corner of the 2x2 sub-grid
                initial_marks.append((i + 1, j + 1))

    # Step 2: Apply transformations
    # First, reflect over the line y=x. A point (r, c) becomes (c, r).
    reflected_marks = [(c, r) for r, c in initial_marks]

    # Second, rotate 90 degrees clockwise. A point (r, c) becomes (c, N-1-r).
    final_marks = [(c, N - 1 - r) for r, c in reflected_marks]

    # Step 3: Create the final grid from the list of marked cells
    final_grid = np.zeros((N, N), dtype=int)
    for r, c in final_marks:
        final_grid[r, c] = 1

    print("The final grid after all transformations is:")
    print(final_grid)
    print("-" * 30)

    # Step 4: Count marked cells in all possible 4x4 sub-grids
    subgrids_with_two_marks = 0
    num_subgrids_to_check = N - subgrid_size + 1
    
    print(f"Checking all {num_subgrids_to_check * num_subgrids_to_check} possible 4x4 sub-grids...")
    
    for R in range(num_subgrids_to_check):
        for C in range(num_subgrids_to_check):
            # Extract the 4x4 sub-grid view
            sub_grid = final_grid[R:R + subgrid_size, C:C + subgrid_size]
            
            # Count the number of marked cells (sum of the elements)
            marked_cells_in_subgrid = np.sum(sub_grid)
            
            # Check if the count is exactly two
            if marked_cells_in_subgrid == 2:
                subgrids_with_two_marks += 1

    # Step 5: Output the final result
    print("-" * 30)
    # The final equation gives the total count.
    print(f"The number of 4x4 sub-grids with exactly 2 marked cells = {subgrids_with_two_marks}")


# Execute the function to find the answer
solve_grid_problem()