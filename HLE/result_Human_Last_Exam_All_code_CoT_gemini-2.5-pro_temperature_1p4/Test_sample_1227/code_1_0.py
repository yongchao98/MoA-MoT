import collections

def solve_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules.
    This solution assumes standard crossword construction conventions, including
    a fixed black border and the rules specified in the problem.
    """
    N = 8
    
    def check_word_length(grid):
        """Checks for any words (white segments) of length 1 or 2."""
        for r in range(N):
            row_str = "".join(map(str, grid[r]))
            for segment in row_str.split('1'):
                if 0 < len(segment) < 3:
                    return False
        for c in range(N):
            col_str = "".join(str(grid[r][c]) for r in range(N))
            for segment in col_str.split('1'):
                if 0 < len(segment) < 3:
                    return False
        return True

    def check_cheaters(grid):
        """
        Checks for "cheater" squares in the inner 6x6 grid. A black square
        is a cheater if it does not separate white squares either horizontally or vertically.
        """
        for r in range(1, N - 1):
            for c in range(1, N - 1):
                if grid[r][c] == 1:
                    is_horiz_splitter = (grid[r][c - 1] == 0 and grid[r][c + 1] == 0)
                    is_vert_splitter = (grid[r - 1][c] == 0 and grid[r + 1][c] == 0)
                    if not is_horiz_splitter and not is_vert_splitter:
                        return False
        return True

    def check_connectivity(grid):
        """Checks if all white squares form a single connected component."""
        white_squares = []
        for r in range(N):
            for c in range(N):
                if grid[r][c] == 0:
                    white_squares.append((r, c))

        if not white_squares:
            return False

        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        count_connected = 1

        while q:
            r, c = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if 0 <= nr < N and 0 <= nc < N and grid[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
                    count_connected += 1
        
        return count_connected == len(white_squares)

    # Define the 18 independent cells in the top half of the inner 6x6 grid
    independent_cells_coords = []
    for r in range(1, N // 2):
        for c in range(1, N - 1):
            independent_cells_coords.append((r, c))
    
    valid_grid_count = 0
    num_independent_cells = len(independent_cells_coords)

    # Iterate through all 2^18 possible patterns
    for i in range(2**num_independent_cells):
        grid = [[1] * N for _ in range(N)]  # Start with an all-black border

        # Populate the grid based on the bits of 'i' (0=white, 1=black)
        for k in range(num_independent_cells):
            if not ((i >> k) & 1):  # If bit is 0, place a white square
                r, c = independent_cells_coords[k]
                grid[r][c] = 0
                grid[N - 1 - r][N - 1 - c] = 0  # Apply rotational symmetry

        # Validate the generated grid
        if check_word_length(grid) and check_cheaters(grid) and check_connectivity(grid):
            valid_grid_count += 1
            
    print(valid_grid_count)

# Execute the solver
solve_crossword_grids()