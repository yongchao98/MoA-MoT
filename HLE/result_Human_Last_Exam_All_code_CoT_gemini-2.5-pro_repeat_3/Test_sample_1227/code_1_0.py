import collections

def count_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules,
    assuming a solid black border around the grid.
    """
    SIZE = 8
    count = 0
    
    # The inner grid is 6x6. Due to rotational symmetry, we only need to decide
    # the colors for the top half of this inner area, which consists of
    # 3 rows of 6 squares, totaling 18 independent cells.
    # The number of combinations to check is 2^18.
    num_combinations = 1 << ((SIZE - 2) * (SIZE - 2) // 2)

    for i in range(num_combinations):
        # Start with a grid that has a solid black border.
        grid = [[1] * SIZE for _ in range(SIZE)]
        
        # Use the bits of 'i' to determine the colors of the 18 independent cells.
        temp_i = i
        # Iterate through the coordinates of the independent cells (top half of inner grid)
        for r_coord in range(1, SIZE // 2):  # Rows 1, 2, 3
            for c_coord in range(1, SIZE - 1): # Columns 1 through 6
                color = temp_i & 1  # 0 for white, 1 for black
                
                # Set the color for the cell and its symmetric counterpart
                grid[r_coord][c_coord] = color
                grid[SIZE - 1 - r_coord][SIZE - 1 - c_coord] = color
                
                temp_i >>= 1
        
        # --- Validation Step ---
        
        # 1. Check for valid word lengths (must be >= 3)
        valid_len = True
        # Check horizontal words
        for r in range(SIZE):
            length = 0
            for c in range(SIZE):
                if grid[r][c] == 0:
                    length += 1
                else:
                    if 0 < length < 3:
                        valid_len = False
                        break
                    length = 0
            if not valid_len or (0 < length < 3):
                valid_len = False
                break
        if not valid_len:
            continue

        # Check vertical words
        for c in range(SIZE):
            length = 0
            for r in range(SIZE):
                if grid[r][c] == 0:
                    length += 1
                else:
                    if 0 < length < 3:
                        valid_len = False
                        break
                    length = 0
            if not valid_len or (0 < length < 3):
                valid_len = False
                break
        if not valid_len:
            continue

        # 2. Check for "cheater" squares (no 2x2 blocks of black squares)
        has_2x2_cheater = False
        for r in range(SIZE - 1):
            for c in range(SIZE - 1):
                if (grid[r][c] == 1 and grid[r + 1][c] == 1 and
                    grid[r][c + 1] == 1 and grid[r + 1][c + 1] == 1):
                    has_2x2_cheater = True
                    break
            if has_2x2_cheater:
                break
        if has_2x2_cheater:
            continue

        # 3. Check for full connectivity of white squares
        white_squares = []
        for r in range(SIZE):
            for c in range(SIZE):
                if grid[r][c] == 0:
                    white_squares.append((r, c))
        
        # If there are no white squares, it's not a valid puzzle.
        if not white_squares:
            continue
        
        # Perform a Breadth-First Search (BFS) to check connectivity.
        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        
        while q:
            r_bfs, c_bfs = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r_bfs + dr, c_bfs + dc
                if (0 <= nr < SIZE and 0 <= nc < SIZE and
                        grid[nr][nc] == 0 and (nr, nc) not in visited):
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        # If the number of visited white squares equals the total, the grid is connected.
        if len(visited) == len(white_squares):
            count += 1
            
    print(count)

if __name__ == '__main__':
    count_crossword_grids()