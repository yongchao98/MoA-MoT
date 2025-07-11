import collections

def solve():
    """
    Solves the grid transformation and counting problem.
    """
    N = 8
    
    # 1. Construct the Initial Grid
    initial_grid = [[0] * N for _ in range(N)]
    num_subgrids_dim = N // 2
    
    print("Step 1: Constructing the initial 8x8 grid.")
    print("The grid is composed of 2x2 sub-grids.")
    print("In the 4x4 grid of sub-grids, the marked cell is at the top-left if (i+j) is even, and bottom-right if (i+j) is odd.")
    
    marked_cells_initial = []
    for i in range(num_subgrids_dim):
        for j in range(num_subgrids_dim):
            # Top-left corner of the 2x2 sub-grid
            base_row, base_col = 2 * i, 2 * j
            
            if (i + j) % 2 == 0:
                # Mark top-left cell
                r, c = base_row, base_col
                initial_grid[r][c] = 1
                marked_cells_initial.append((r,c))
            else:
                # Mark bottom-right cell
                r, c = base_row + 1, base_col + 1
                initial_grid[r][c] = 1
                marked_cells_initial.append((r,c))

    # 2. Apply Transformations
    # 2a. Reflect over y=x (Transpose)
    print("\nStep 2: Applying transformations.")
    print("Reflecting the grid over the line y=x (transposing the matrix).")
    reflected_grid = [[0] * N for _ in range(N)]
    for r in range(N):
        for c in range(N):
            reflected_grid[c][r] = initial_grid[r][c]

    # 2b. Rotate 90 degrees clockwise
    print("Rotating the resulting grid 90 degrees clockwise.")
    final_grid = [[0] * N for _ in range(N)]
    for r in range(N):
        for c in range(N):
            final_grid[c][N - 1 - r] = reflected_grid[r][c]
            
    # 3. Count 4x4 Sub-grids
    print("\nStep 3: Counting 4x4 sub-grids with exactly two marked cells.")
    subgrid_size = 4
    count_with_two_marked = 0
    
    # There are (N - subgrid_size + 1) * (N - subgrid_size + 1) possible sub-grids
    num_subgrids_to_check = N - subgrid_size + 1
    
    subgrid_sums = collections.defaultdict(int)

    for r_start in range(num_subgrids_to_check):
        for c_start in range(num_subgrids_to_check):
            marked_count = 0
            for r in range(subgrid_size):
                for c in range(subgrid_size):
                    if final_grid[r_start + r][c_start + c] == 1:
                        marked_count += 1
            
            subgrid_sums[marked_count] += 1
            if marked_count == 2:
                count_with_two_marked += 1

    print("\nAnalysis of all 25 possible 4x4 sub-grids:")
    if not subgrid_sums:
        print("No sub-grids found.")
    else:
        for num_marked, count in sorted(subgrid_sums.items()):
            print(f"- Found {count} sub-grid(s) with {num_marked} marked cells.")

    print("\nFinal Calculation:")
    total_subgrids_with_2_marked = 0
    equation_parts = []
    for num_marked, count in sorted(subgrid_sums.items()):
        if num_marked == 2:
             total_subgrids_with_2_marked += count
             equation_parts.append(str(count))
    
    if not equation_parts:
        equation_parts.append('0')

    # The problem asks for the maximum number. Since the configuration is fixed,
    # this is just the total number of such sub-grids that exist.
    print(f"The number of 4x4 sub-grids with exactly two marked cells is the sum of all such grids found.")
    print(f"Total = {' + '.join(equation_parts)} = {total_subgrids_with_2_marked}")
    
    print("\nThe final answer is the total count.")
    print(f"{total_subgrids_with_2_marked}")

solve()
<<<0>>>