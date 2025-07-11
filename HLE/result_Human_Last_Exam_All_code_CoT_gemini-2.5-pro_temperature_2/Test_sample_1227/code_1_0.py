import sys

# The recursive search can go deep, so we increase the recursion limit.
sys.setrecursionlimit(2000)

SIZE = 8
# A global variable for the grid makes the recursive function simpler.
# It stores the grid state: 0 for white, 1 for black, -1 for undecided.
grid = [[-1] * SIZE for _ in range(SIZE)]
# The list of the 32 cells that, due to symmetry, determine the whole grid.
# We will iterate through these cells to build the grid.
cells_to_set = []
for r in range(SIZE // 2):
    for c in range(SIZE):
        cells_to_set.append((r, c))

# --- Constraint Checking Helper Functions ---

def check_line_word_length(line):
    """Checks a single row or column for words shorter than 3 letters."""
    # We pad the line with black squares (1) to easily check words at the edges.
    padded_line = [1] + line + [1]
    # Look for a pattern of a single white square between two black squares (BWB).
    for i in range(len(padded_line) - 2):
        if padded_line[i] == 1 and padded_line[i+1] == 0 and padded_line[i+2] == 1:
            return False
    # Look for a pattern of two white squares between two black squares (BWWB).
    for i in range(len(padded_line) - 3):
        if padded_line[i] == 1 and padded_line[i+1] == 0 and padded_line[i+2] == 0 and padded_line[i+3] == 1:
            return False
    return True

def check_all_columns():
    """Applies the word length check to all columns of the grid."""
    for c in range(SIZE):
        col = [grid[r][c] for r in range(SIZE)]
        if not check_line_word_length(col):
            return False
    return True

def is_connected():
    """Checks if all white squares on the grid form a single connected component."""
    white_squares = []
    for r in range(SIZE):
        for c in range(SIZE):
            if grid[r][c] == 0:
                white_squares.append((r, c))

    # If there are no white squares, it's considered connected.
    # However, such grids will be eliminated by the "no cheaters" rule.
    if not white_squares:
        return True

    # Perform a breadth-first search (BFS) starting from the first white square found.
    total_white = len(white_squares)
    q = [white_squares[0]]
    visited = {white_squares[0]}
    count = 0
    while q:
        r, c = q.pop(0)
        count += 1
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if (0 <= nr < SIZE and 0 <= nc < SIZE and
                grid[nr][nc] == 0 and (nr, nc) not in visited):
                visited.add((nr, nc))
                q.append((nr, nc))
    
    # If the number of visited squares equals the total number, the grid is connected.
    return count == total_white

def has_no_cheaters():
    """Checks for "cheater" squares (black squares not adjacent to any white square)."""
    for r in range(SIZE):
        for c in range(SIZE):
            if grid[r][c] == 1:  # It's a black square
                is_cheater = True
                # Check all four neighbors.
                for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                    nr, nc = r + dr, c + dc
                    # If a neighbor is inside the grid and is white, it's not a cheater.
                    if (0 <= nr < SIZE and 0 <= nc < SIZE and grid[nr][nc] == 0):
                        is_cheater = False
                        break
                if is_cheater:
                    return False
    return True

# --- The Backtracking Solver ---

# This global variable will hold the final count of valid grids.
final_grid_count = 0

def backtrack_solver(k):
    """
    Recursively fills the grid, pruning branches that violate constraints.
    'k' is the index into the `cells_to_set` list, tracking our progress.
    """
    global final_grid_count
    
    # Base Case: If we have decided all 32 independent cells, the grid is complete.
    if k == len(cells_to_set):
        # By this point, the grid is already confirmed to be symmetric and
        # to have no short words. We just need to check the final two constraints.
        if is_connected() and has_no_cheaters():
            final_grid_count += 1
        return

    r, c = cells_to_set[k]
    r_sym, c_sym = SIZE - 1 - r, SIZE - 1 - c

    # --- Recursive Step: Try both white and black for the current cell ---

    # Option 1: Place a WHITE square (0)
    grid[r][c] = 0
    grid[r_sym][c_sym] = 0
    
    proceed = True
    # Pruning Check: If a row is complete, check its word lengths immediately.
    if c == SIZE - 1: 
        if not (check_line_word_length(grid[r]) and check_line_word_length(grid[r_sym])):
            proceed = False
        # If we just finished the last independent row, the grid is full. Check columns.
        if proceed and r == (SIZE // 2) - 1:
            if not check_all_columns():
                proceed = False
    if proceed:
        backtrack_solver(k + 1)

    # Option 2: Place a BLACK square (1)
    grid[r][c] = 1
    grid[r_sym][c_sym] = 1
    
    proceed = True
    if c == SIZE - 1: # Row is complete, check word lengths.
        if not (check_line_word_length(grid[r]) and check_line_word_length(grid[r_sym])):
            proceed = False
        if proceed and r == (SIZE // 2) - 1: # Grid is full, check columns.
            if not check_all_columns():
                proceed = False
    if proceed:
        backtrack_solver(k + 1)
        
# --- Main Execution ---
# Kick off the recursive search from the first independent cell (k=0).
# The calculation might take several seconds to a few minutes to complete.
backtrack_solver(0)

# The user prompt mentions "output each number in the final equation!".
# This doesn't directly apply here, as we are solving for a single total count.
# So we will just print the final number, which is the answer to the user's question.
print(final_grid_count)