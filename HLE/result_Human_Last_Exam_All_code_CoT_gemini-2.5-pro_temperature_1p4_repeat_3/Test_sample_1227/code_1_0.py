import collections

# This script calculates the number of valid 8x8 crossword grids based on a set of rules.
# It uses a backtracking search algorithm to explore all possible grid patterns that respect
# 180-degree rotational symmetry and a pre-filled black border.

# Global state variables for the search
SIZE = 8
BLACK = 1
WHITE = 0
grid = [[WHITE] * SIZE for _ in range(SIZE)]
solution_count = 0

def check_no_2x2_black():
    """Checks for any 2x2 blocks of black squares in the grid."""
    for r in range(SIZE - 1):
        for c in range(SIZE - 1):
            if (grid[r][c] == BLACK and
                grid[r + 1][c] == BLACK and
                grid[r][c + 1] == BLACK and
                grid[r + 1][c + 1] == BLACK):
                return False
    return True

def check_word_lengths():
    """Ensures all horizontal and vertical runs of white squares (words) are of length 3 or more."""
    # Check rows
    for r in range(SIZE):
        length = 0
        for c in range(SIZE):
            if grid[r][c] == WHITE:
                length += 1
            else:
                if 0 < length < 3:
                    return False
                length = 0
        if 0 < length < 3:
            return False
    # Check columns
    for c in range(SIZE):
        length = 0
        for r in range(SIZE):
            if grid[r][c] == WHITE:
                length += 1
            else:
                if 0 < length < 3:
                    return False
                length = 0
        if 0 < length < 3:
            return False
    return True

def check_connectivity():
    """Verifies that all white squares in the grid form a single connected component using BFS."""
    white_squares = []
    for r in range(SIZE):
        for c in range(SIZE):
            if grid[r][c] == WHITE:
                white_squares.append((r, c))

    if not white_squares:
        # An all-black grid will be invalid due to 2x2 blocks, but is technically "connected".
        return True

    # BFS to find all connected white squares from the first one found
    q = collections.deque([white_squares[0]])
    visited = {white_squares[0]}
    
    count = 0
    while q:
        r, c = q.popleft()
        count += 1
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if (0 <= nr < SIZE and 0 <= nc < SIZE and
                grid[nr][nc] == WHITE and (nr, nc) not in visited):
                visited.add((nr, nc))
                q.append((nr, nc))

    return count == len(white_squares)

def is_valid():
    """Checks if a completed grid satisfies all rules."""
    if not check_no_2x2_black():
        return False
    if not check_word_lengths():
        return False
    if not check_connectivity():
        return False
    return True

def search(k, independent_cells):
    """
    Recursively explores grid configurations.
    
    k: The index of the current independent cell to decide.
    independent_cells: A list of coordinates for the cells that need to be decided.
    """
    global solution_count
    
    if k == len(independent_cells):
        if is_valid():
            solution_count += 1
        return

    r, c = independent_cells[k]
    r_s, c_s = SIZE - 1 - r, SIZE - 1 - c

    # Branch 1: Set the symmetric pair of cells to WHITE
    grid[r][c] = WHITE
    grid[r_s][c_s] = WHITE
    search(k + 1, independent_cells)

    # Branch 2: Set the symmetric pair of cells to BLACK
    grid[r][c] = BLACK
    grid[r_s][c_s] = BLACK
    search(k + 1, independent_cells)

def solve_puzzle():
    """Initializes the grid and starts the backtracking search."""
    
    # Assumption: The outermost border is black.
    for i in range(SIZE):
        grid[0][i] = BLACK
        grid[SIZE - 1][i] = BLACK
        grid[i][0] = BLACK
        grid[i][SIZE - 1] = BLACK

    # Due to symmetry, we only need to decide the state of the top half of the inner 6x6 grid.
    # These are 18 cells, which will determine the entire grid pattern.
    independent_cells = []
    for r in range(1, SIZE // 2):      # Rows 1, 2, 3
        for c in range(1, SIZE - 1):   # Columns 1, 2, 3, 4, 5, 6
            independent_cells.append((r,c))
            
    # Start the search from the first independent cell.
    search(0, independent_cells)
    
    print(solution_count)

if __name__ == '__main__':
    solve_puzzle()