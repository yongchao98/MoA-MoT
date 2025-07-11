import collections

def solve_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules.
    This script iterates through all possible symmetric patterns for the inner
    6x6 grid and validates each one against the specified constraints.
    """
    # The final count of the equation we are solving is the number of valid grids.
    # The equation is essentially: Count(all grids) where Grid satisfies all rules.
    
    # Grid parameters
    N = 8
    valid_grid_count = 0

    # Because of 180-degree symmetry and an assumed black border, we only need
    # to define the state of the top half of the inner 6x6 grid. This consists
    # of 18 cells, which we define here.
    independent_cells = []
    for r in range(1, 4):  # Rows 1, 2, 3 of the 8x8 grid
        for c in range(1, 7):  # Columns 1 through 6
            independent_cells.append((r, c))
    
    num_independent = len(independent_cells) # This will be 18

    # Iterate through all 2^18 = 262,144 possible patterns for the independent cells
    for i in range(1 << num_independent):
        # 1. Build the grid for the current pattern `i`
        grid = [[-1] * N for _ in range(N)]
        
        # Assume black borders
        for k in range(N):
            grid[0][k] = grid[N - 1][k] = 1
            grid[k][0] = grid[k][N - 1] = 1

        # Use the bits of `i` to fill the independent cells and their symmetric pairs
        temp_i = i
        for r_cell, c_cell in independent_cells:
            color = temp_i & 1  # 0 for white, 1 for black
            grid[r_cell][c_cell] = color
            grid[N - 1 - r_cell][N - 1 - c_cell] = color
            temp_i >>= 1
        
        # 2. Validate the generated grid
        
        # Check for 2x2 black squares ("cheaters")
        has_2x2_black = False
        for r in range(N - 1):
            for c in range(N - 1):
                if (grid[r][c] == 1 and grid[r + 1][c] == 1 and
                    grid[r][c + 1] == 1 and grid[r + 1][c + 1] == 1):
                    has_2x2_black = True
                    break
            if has_2x2_black:
                break
        if has_2x2_black:
            continue # Invalid grid, try next pattern

        # Check word lengths (min 3)
        invalid_length = False
        for k in range(N):
            # Check row k
            row = grid[k]
            length = 0
            padded_row = [1] + row + [1]
            for cell in padded_row:
                if cell == 0:
                    length += 1
                else:
                    if 0 < length < 3:
                        invalid_length = True
                        break
                    length = 0
            if invalid_length: break

            # Check column k
            col = [grid[r][k] for r in range(N)]
            length = 0
            padded_col = [1] + col + [1]
            for cell in padded_col:
                if cell == 0:
                    length += 1
                else:
                    if 0 < length < 3:
                        invalid_length = True
                        break
                    length = 0
            if invalid_length: break
        
        if invalid_length:
            continue # Invalid grid, try next pattern

        # Check for full interconnectivity of white squares
        white_squares = []
        for r in range(N):
            for c in range(N):
                if grid[r][c] == 0:
                    white_squares.append((r, c))
        
        total_white_squares = len(white_squares)
        
        if total_white_squares == 0:
            continue # A grid with no white squares isn't a valid puzzle

        # Perform a BFS to count connected white squares
        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        
        while q:
            r_q, c_q = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r_q + dr, c_q + dc
                if (0 <= nr < N and 0 <= nc < N and 
                    grid[nr][nc] == 0 and (nr, nc) not in visited):
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        # If all white squares were visited, the grid is fully connected
        if len(visited) == total_white_squares:
            valid_grid_count += 1
    
    # Final Output: Print the numbers involved in finding the final answer.
    print(f"Number of independent cells to check: {num_independent}")
    print(f"Total symmetric patterns to test: 2^{num_independent} = {1 << num_independent}")
    print("\n------------------------------------------------------")
    print("Final Answer:")
    print(f"The number of possible 8x8 crossword grids is {valid_grid_count}")

# Execute the solver
solve_crossword_grids()