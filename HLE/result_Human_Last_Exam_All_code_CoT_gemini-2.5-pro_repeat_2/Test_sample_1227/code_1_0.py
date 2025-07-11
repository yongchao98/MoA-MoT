import collections

def count_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules.
    
    The rules are:
    1. 8x8 grid size.
    2. 180-degree rotational symmetry.
    3. Minimum word length of 3 (implies black border).
    4. All white squares are fully connected.
    5. No "cheater" squares (interpreted as no 2x2 black squares).
    """
    N = 8
    valid_grid_count = 0

    # We only need to choose the color for the 18 unique squares in the top half
    # of the inner 6x6 grid. The rest of the grid is determined by symmetry
    # and the black border rule.
    cells_to_iterate = []
    for r in range(1, N // 2):      # Rows 1, 2, 3
        for c in range(1, N - 1):  # Columns 1 to 6
            cells_to_iterate.append((r, c))
    
    num_deciding_cells = len(cells_to_iterate)

    # Iterate through all 2^18 possible patterns for the inner grid
    for i in range(1 << num_deciding_cells):
        # 1. Construct the grid for the current pattern
        grid = [[0] * N for _ in range(N)] # 0 for white, 1 for black

        # Rule: The border is always black
        for k in range(N):
            grid[0][k] = 1
            grid[N - 1][k] = 1
            grid[k][0] = 1
            grid[k][N - 1] = 1

        # Fill the inner grid based on the bits of i, applying symmetry
        temp_i = i
        for r, c in cells_to_iterate:
            color = temp_i & 1
            grid[r][c] = color
            grid[N - 1 - r][N - 1 - c] = color  # Apply 180-degree symmetry
            temp_i >>= 1

        # 2. Validate the generated grid against the rules
        is_valid = True

        # Rule: No 2x2 blocks of black squares
        for r in range(N - 1):
            for c in range(N - 1):
                if (grid[r][c] == 1 and grid[r + 1][c] == 1 and
                        grid[r][c + 1] == 1 and grid[r + 1][c + 1] == 1):
                    is_valid = False
                    break
            if not is_valid: break
        if not is_valid: continue

        # Rule: Minimum word length of 3 (check for invalid short words)
        for r in range(N): # Check all rows
            for c in range(N - 2): # Check for 'BWB' pattern (1-letter word)
                if grid[r][c] == 1 and grid[r][c + 1] == 0 and grid[r][c + 2] == 1:
                    is_valid = False; break
            if not is_valid: break
            for c in range(N - 3): # Check for 'BWWB' pattern (2-letter word)
                if grid[r][c] == 1 and grid[r][c+1] == 0 and grid[r][c+2] == 0 and grid[r][c+3] == 1:
                    is_valid = False; break
            if not is_valid: break
        if not is_valid: continue

        for c in range(N): # Check all columns
            for r in range(N - 2): # Check for 'BWB' pattern
                if grid[r][c] == 1 and grid[r + 1][c] == 0 and grid[r + 2][c] == 1:
                    is_valid = False; break
            if not is_valid: break
            for r in range(N - 3): # Check for 'BWWB' pattern
                if grid[r][c] == 1 and grid[r+1][c] == 0 and grid[r+2][c] == 0 and grid[r+3][c] == 1:
                    is_valid = False; break
            if not is_valid: break
        if not is_valid: continue

        # Rule: All white squares must be fully connected
        white_squares = []
        for r in range(N):
            for c in range(N):
                if grid[r][c] == 0:
                    white_squares.append((r, c))

        if not white_squares: # A valid puzzle must have words.
            continue
        
        # Use Breadth-First Search (BFS) to check connectivity
        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        while q:
            row, col = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = row + dr, col + dc
                if 0 <= nr < N and 0 <= nc < N and grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        if len(visited) != len(white_squares):
            continue

        # 3. If all rules are met, this is a valid grid
        valid_grid_count += 1
        
    print(valid_grid_count)

if __name__ == '__main__':
    count_crossword_grids()