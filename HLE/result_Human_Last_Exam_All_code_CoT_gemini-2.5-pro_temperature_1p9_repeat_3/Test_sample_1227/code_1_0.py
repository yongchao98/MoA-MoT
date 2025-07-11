import sys

def solve_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules.

    The solution assumes a standard convention that the grid is framed by black
    squares, reducing the problem to filling the inner 6x6 area while respecting
    the grid's 180-degree rotational symmetry.
    """
    # Set a higher recursion limit, just in case, though 18 levels is shallow.
    sys.setrecursionlimit(2000)

    N = 8
    SUB_N = 6
    
    # Initialize an 8x8 grid with a black frame. 1=Black, 0=White.
    grid = [[1] * N for _ in range(N)]

    # Define the 18 independent squares within the inner 6x6 subgrid.
    # We select the top half of the subgrid's rows.
    # Rows 1, 2, 3 of the 8x8 grid.
    independent_squares = []
    for r_sub in range(SUB_N // 2):      # r_sub is 0, 1, 2
        for c_sub in range(SUB_N):       # c_sub is 0, 1, 2, 3, 4, 5
            # Map subgrid coordinates to the main 8x8 grid coordinates
            r, c = r_sub + 1, c_sub + 1
            independent_squares.append((r, c))
    
    # Find the number of valid grids using recursion
    total_valid_grids = find_grids_recursive(0, independent_squares, grid, N)
    
    print(f"The number of possible 8x8 crossword grids is: {total_valid_grids}")

def find_grids_recursive(k, independent_squares, grid, N):
    """
    Recursively explores all possible grid colorings and returns the count of valid ones.
    """
    # Base case: All independent squares have been colored.
    if k == len(independent_squares):
        # Validate the complete grid against the rules.
        if is_valid_grid(grid, N):
            return 1
        else:
            return 0

    # Get the coordinates for the current independent square and its symmetric partner.
    r, c = independent_squares[k]
    r_sym, c_sym = N - 1 - r, N - 1 - c

    total_count = 0
    
    # --- Try coloring the pair of squares WHITE (0) ---
    grid[r][c] = 0
    grid[r_sym][c_sym] = 0
    total_count += find_grids_recursive(k + 1, independent_squares, grid, N)
    
    # --- Try coloring the pair of squares BLACK (1) ---
    grid[r][c] = 1
    grid[r_sym][c_sym] = 1
    total_count += find_grids_recursive(k + 1, independent_squares, grid, N)
    
    return total_count

def is_valid_grid(grid, N):
    """Checks if a grid satisfies all the required crossword puzzle rules."""
    return (check_word_lengths(grid, N) and
            check_connectivity(grid, N) and
            check_no_cheaters(grid, N))

def check_word_lengths(grid, N):
    """Ensures all words (runs of white squares) are at least 3 letters long."""
    # Check rows
    for r in range(N):
        length = 0
        for c in range(N):
            if grid[r][c] == 0:  # White square
                length += 1
            else:  # Black square
                if 1 <= length <= 2:
                    return False
                length = 0
        if 1 <= length <= 2: # Check run at the end of the row
            return False
            
    # Check columns
    for c in range(N):
        length = 0
        for r in range(N):
            if grid[r][c] == 0:
                length += 1
            else:
                if 1 <= length <= 2:
                    return False
                length = 0
        if 1 <= length <= 2: # Check run at the end of the column
            return False
            
    return True

def check_connectivity(grid, N):
    """Ensures all white squares are connected into a single component."""
    white_squares = []
    for r in range(N):
        for c in range(N):
            if grid[r][c] == 0:
                white_squares.append((r, c))
    
    if not white_squares:
        return False  # An all-black grid is not a valid puzzle.

    # BFS traversal to find all reachable white squares from the first one.
    q = [white_squares[0]]
    visited = {white_squares[0]}
    head = 0
    while head < len(q):
        r, c = q[head]
        head += 1
        
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if (0 <= nr < N and 0 <= nc < N and 
                    grid[nr][nc] == 0 and (nr, nc) not in visited):
                visited.add((nr, nc))
                q.append((nr, nc))
                
    return len(visited) == len(white_squares)

def check_no_cheaters(grid, N):
    """Ensures there are no 2x2 blocks of black squares."""
    for r in range(N - 1):
        for c in range(N - 1):
            if (grid[r][c] == 1 and grid[r+1][c] == 1 and
                    grid[r][c+1] == 1 and grid[r+1][c+1] == 1):
                return False
    return True

if __name__ == '__main__':
    solve_crossword_grids()