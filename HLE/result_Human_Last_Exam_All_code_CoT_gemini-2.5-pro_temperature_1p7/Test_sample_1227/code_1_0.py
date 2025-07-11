import collections

def solve_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a specific set of rules.
    """
    N = 8
    
    # We use a mutable list to hold the count so it can be modified by the recursive function.
    count = [0]
    
    # Grid: -1: undecided, 0: white, 1: black
    grid = [[-1] * N for _ in range(N)]
    
    # Independent cells in the top half of the grid that we will decide.
    independent_cells = []
    for r in range(N // 2):
        for c in range(N):
            independent_cells.append((r, c))

    def check_cheaters(current_grid):
        """
        Checks for "cheater" squares (white squares with >= 3 black neighbors).
        Boundaries count as black neighbors. Valid for partial grids.
        """
        for r in range(N):
            for c in range(N):
                if current_grid[r][c] == 0:  # If it's a white square
                    black_neighbors = 0
                    for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                        nr, nc = r + dr, c + dc
                        if not (0 <= nr < N and 0 <= nc < N) or current_grid[nr][nc] == 1:
                            black_neighbors += 1
                    if black_neighbors >= 3:
                        return False
        return True

    def check_short_words(current_grid):
        """
        Checks for completed words with length < 3.
        Valid for partial grids as it only checks definitively terminated words.
        """
        # Horizontal check
        for r in range(N):
            for c in range(N):
                if current_grid[r][c] == 0 and (c == 0 or current_grid[r][c-1] == 1):
                    word_len = 0
                    is_terminated = False
                    for k in range(c, N):
                        if current_grid[r][k] == 0:
                            word_len += 1
                        elif current_grid[r][k] == 1:
                            is_terminated = True
                            break
                        else:  # Undecided square, word is not terminated yet.
                            break
                    else: # Reached end of row, word is terminated.
                        is_terminated = True
                    
                    if is_terminated and word_len < 3:
                        return False

        # Vertical check
        for c in range(N):
            for r in range(N):
                if current_grid[r][c] == 0 and (r == 0 or current_grid[r-1][c] == 1):
                    word_len = 0
                    is_terminated = False
                    for k in range(r, N):
                        if current_grid[k][c] == 0:
                            word_len += 1
                        elif current_grid[k][c] == 1:
                            is_terminated = True
                            break
                        else: # Undecided square
                            break
                    else: # Reached end of column
                        is_terminated = True

                    if is_terminated and word_len < 3:
                        return False
        return True

    def check_connectivity(final_grid):
        """
        Checks if all white squares are connected. Runs only on a full grid.
        """
        white_squares = []
        for r in range(N):
            for c in range(N):
                if final_grid[r][c] == 0:
                    white_squares.append((r,c))

        if not white_squares: # All-black grid is valid.
            return True

        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        
        while q:
            r, c = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if (0 <= nr < N and 0 <= nc < N and
                        final_grid[nr][nc] == 0 and (nr, nc) not in visited):
                    visited.add((nr, nc))
                    q.append((nr, nc))

        return len(visited) == len(white_squares)

    def solve_recursive(cell_index):
        """
        Recursively explores grid configurations with backtracking.
        """
        if cell_index == len(independent_cells):
            if check_connectivity(grid):
                count[0] += 1
            return

        r, c = independent_cells[cell_index]
        sym_r, sym_c = N - 1 - r, N - 1 - c

        # Try placing a WHITE square
        grid[r][c] = 0
        grid[sym_r][sym_c] = 0
        if check_cheaters(grid):
            solve_recursive(cell_index + 1)

        # Try placing a BLACK square
        grid[r][c] = 1
        grid[sym_r][sym_c] = 1
        if check_cheaters(grid) and check_short_words(grid):
            solve_recursive(cell_index + 1)
        
        # Backtrack (restore grid state for parent call)
        grid[r][c] = -1
        grid[sym_r][sym_c] = -1

    solve_recursive(0)
    print(count[0])

if __name__ == '__main__':
    solve_crossword_grids()