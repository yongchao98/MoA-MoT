import numpy as np

def solve_grid_problem():
    """
    This function solves the grid transformation problem by:
    1. Creating the initial 8x8 grid with marked cells based on the checkerboard pattern.
    2. Applying the transformations (reflection and rotation), which simplifies to a horizontal flip.
    3. Iterating through all 25 possible 4x4 sub-grids in the final grid.
    4. Counting the number of marked cells in each sub-grid.
    5. Reporting the number of sub-grids that have exactly two marked cells.
    """
    N = 8
    
    # 1. Create and populate the initial grid
    # We choose one of the two possible checkerboard patterns. The other pattern gives an
    # isomorphic result, so the "maximum" is not affected by this choice.
    initial_grid = np.zeros((N, N), dtype=int)
    for i in range(N // 2):  # 2x2 sub-grid row index (0-3)
        for j in range(N // 2):  # 2x2 sub-grid column index (0-3)
            base_r, base_c = 2 * i, 2 * j
            # Checkerboard pattern for which corner to mark
            if (i + j) % 2 == 0:
                initial_grid[base_r, base_c] = 1  # Top-left corner
            else:
                initial_grid[base_r + 1, base_c + 1] = 1  # Bottom-right corner

    # 2. Apply transformations
    # Reflection over y=x (transpose) and then 90-degree clockwise rotation
    # is equivalent to a horizontal flip (fliplr).
    final_grid = np.fliplr(initial_grid)

    # 3. & 4. Iterate through 4x4 sub-grids and find their sums
    subgrid_size = 4
    count_with_two_marked = 0
    
    print("Number of marked cells in each 4x4 sub-grid:")
    equation_parts = []

    num_subgrid_rows = N - subgrid_size + 1
    num_subgrid_cols = N - subgrid_size + 1

    for r in range(num_subgrid_rows):
        for c in range(num_subgrid_cols):
            subgrid = final_grid[r:r + subgrid_size, c:c + subgrid_size]
            current_sum = np.sum(subgrid)
            
            # This part satisfies the "output each number in the final equation" requirement.
            # We'll build a list of all sums to show the components.
            if current_sum == 2:
                equation_parts.append("1") # Add 1 to the sum for each subgrid that matches
            else:
                equation_parts.append("0") # Add 0 if it doesn't match

            print(f"Sub-grid at ({r},{c}): sum = {current_sum}")
            
            if current_sum == 2:
                count_with_two_marked += 1
    
    # 5. Output the final result
    # We display the "equation" as a sum of 1s and 0s
    print("\nTo find the total number of sub-grids with exactly two marked cells, we sum '1' for each match and '0' otherwise:")
    final_equation = " + ".join(equation_parts)
    print(f"{final_equation} = {count_with_two_marked}")
    
    print("\nThe maximum number of 4x4 sub-grids with exactly two marked cells is:")
    print(count_with_two_marked)

solve_grid_problem()