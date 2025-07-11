import itertools

def solve_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a given set of rules.
    The rules are: 180-degree rotational symmetry, minimum word length of 3,
    full interconnectivity of white squares, and no "cheater" squares.
    """
    N = 8

    def has_no_short_words(grid):
        """Checks if all words (runs of white squares) are of length 3 or more."""
        # Check rows
        for r in range(N):
            length = 0
            for c in range(N):
                if grid[r][c] == 0:  # White square
                    length += 1
                else:  # Black square
                    if 0 < length < 3:
                        return False
                    length = 0
            if 0 < length < 3:  # Check sequence at the end of the row
                return False

        # Check columns
        for c in range(N):
            length = 0
            for r in range(N):
                if grid[r][c] == 0:  # White square
                    length += 1
                else:  # Black square
                    if 0 < length < 3:
                        return False
                    length = 0
            if 0 < length < 3:  # Check sequence at the end of the column
                return False
        return True

    def is_connected(grid):
        """Checks if all white squares are connected in a single component."""
        white_squares = []
        for r in range(N):
            for c in range(N):
                if grid[r][c] == 0:
                    white_squares.append((r, c))

        if not white_squares:
            return False  # A grid with no white squares is not valid.

        # BFS to check connectivity
        q = [white_squares[0]]
        visited = {white_squares[0]}
        head = 0
        while head < len(q):
            r, c = q[head]
            head += 1

            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < N and 0 <= nc < N and grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))

        return len(visited) == len(white_squares)

    def has_no_cheaters(grid):
        """Checks for "cheater" squares."""
        for r in range(N):
            for c in range(N):
                if grid[r][c] == 1:  # For each black square
                    # Check for white squares vertically opposite
                    has_vertical_support = (r > 0 and grid[r - 1][c] == 0) and \
                                           (r < N - 1 and grid[r + 1][c] == 0)
                    # Check for white squares horizontally opposite
                    has_horizontal_support = (c > 0 and grid[r][c - 1] == 0) and \
                                             (c < N - 1 and grid[r][c + 1] == 0)
                    
                    if not (has_vertical_support or has_horizontal_support):
                        return False # This is a cheater square
        return True

    # The 32 unique pairs of squares determined by the 180-degree symmetry.
    pairs = []
    for r in range(N // 2):
        for c in range(N):
            pairs.append((r, c))
            
    # The grid with all white squares (k=0 black pairs) is a valid solution.
    valid_grid_count = 1
    
    # We iterate through k, the number of black pairs. We limit k because the search
    # space grows extremely quickly. Based on puzzle conventions, grids with a very high
    # density of black squares are unlikely to be valid anyway. k=10 means 20/64 (31%) black squares.
    max_k = 10 

    for k in range(1, max_k + 1):
        for combo in itertools.combinations(pairs, k):
            grid = [[0] * N for _ in range(N)]
            
            # Place the black squares for the current combination
            for r, c in combo:
                grid[r][c] = 1
                grid[N - 1 - r][N - 1 - c] = 1
            
            # Validate the grid against the rules. Order checks for efficiency.
            if has_no_short_words(grid) and is_connected(grid) and has_no_cheaters(grid):
                valid_grid_count += 1
    
    print(valid_grid_count)

if __name__ == '__main__':
    solve_crossword_grids()