import collections

def solve_crossword_grids():
    """
    This script calculates the number of valid 8x8 crossword grids
    based on a set of rules.

    The plan is to iterate through all possible symmetric patterns for the inner
    6x6 grid and validate each one against the given criteria.
    """
    size = 8
    count = 0
    inner_size = size - 2 # 6 for an 8x8 grid

    # Due to 180-degree symmetry, we only need to decide the color for the
    # top half of the inner grid. For a 6x6 inner grid, this is 18 cells.
    num_independent_cells = (inner_size * inner_size) // 2

    # Pre-calculate the coordinates of these independent cells.
    # These are cells in the top half of the inner grid (rows 1-3, cols 1-6).
    independent_cells_coords = []
    # In a 6x6 inner grid (indices 1 to 6), the top half rows are 1, 2, 3.
    for r in range(1, 1 + inner_size // 2):
        for c in range(1, 1 + inner_size):
            independent_cells_coords.append((r, c))

    # Helper function to check if a square is white (1), treating grid boundaries as black (0).
    def is_white(r, c, grid, grid_size):
        if 0 <= r < grid_size and 0 <= c < grid_size:
            return grid[r][c] == 1
        return False

    def validate(grid, grid_size):
        """Checks if a complete grid is a valid crossword pattern."""
        total_white = sum(row.count(1) for row in grid)
        if total_white == 0:
            return False

        # Rule 1: Minimum word length of 3
        for r in range(grid_size):
            for c in range(grid_size):
                # Check for the start of a horizontal word
                if is_white(r, c, grid, grid_size) and not is_white(r, c - 1, grid, grid_size):
                    length = 0
                    while is_white(r, c + length, grid, grid_size):
                        length += 1
                    if length < 3:
                        return False
                # Check for the start of a vertical word
                if is_white(r, c, grid, grid_size) and not is_white(r - 1, c, grid, grid_size):
                    length = 0
                    while is_white(r + length, c, grid, grid_size):
                        length += 1
                    if length < 3:
                        return False

        # Rule 2: No "cheater" squares
        for r in range(grid_size):
            for c in range(grid_size):
                if grid[r][c] == 0:  # If black square
                    bridges_h = is_white(r, c - 1, grid, grid_size) and is_white(r, c + 1, grid, grid_size)
                    bridges_v = is_white(r - 1, c, grid, grid_size) and is_white(r + 1, c, grid, grid_size)
                    if not bridges_h and not bridges_v:
                        return False
        
        # Rule 3: Full interconnectivity of white squares
        q = collections.deque()
        visited = set()
        first_white_found = False
        for r in range(grid_size):
            for c in range(grid_size):
                if grid[r][c] == 1:
                    q.append((r, c))
                    visited.add((r, c))
                    first_white_found = True
                    break
            if first_white_found:
                break
        
        while q:
            r, c = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if is_white(nr, nc, grid, grid_size) and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        if len(visited) != total_white:
            return False

        return True

    # Iterate through all 2^18 possibilities
    for i in range(2 ** num_independent_cells):
        grid = [[0 for _ in range(size)] for _ in range(size)]
        temp_i = i
        for r_cell, c_cell in independent_cells_coords:
            if (temp_i & 1) == 1:
                grid[r_cell][c_cell] = 1
                r_sym = size - 1 - r_cell
                c_sym = size - 1 - c_cell
                grid[r_sym][c_sym] = 1
            temp_i >>= 1
        
        if validate(grid, size):
            count += 1
            
    print(f"The number of possible 8x8 crossword grids is: {count}")

solve_crossword_grids()