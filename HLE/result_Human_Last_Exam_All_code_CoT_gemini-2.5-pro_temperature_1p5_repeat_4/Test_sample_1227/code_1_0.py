import sys

# Increase recursion limit for deep search on larger grids.
sys.setrecursionlimit(2000)

SIZE = 8
BLACK, WHITE, UNDECIDED = 1, 0, -1

# Global variables to store the grid state and the final count.
grid = [[UNDECIDED] * SIZE for _ in range(SIZE)]
solution_count = 0
# Pre-calculate the cells to iterate over to enforce symmetry.
cells_to_fill = []
for r in range(SIZE):
    for c in range(SIZE):
        # Only decide for the top half of cells.
        if r * SIZE + c < (SIZE * SIZE) / 2:
            cells_to_fill.append((r, c))

def check_connectivity():
    """Checks if all white squares are connected using BFS."""
    white_squares = []
    for r in range(SIZE):
        for c in range(SIZE):
            if grid[r][c] == WHITE:
                white_squares.append((r, c))

    if not white_squares:
        return True # An all-black grid is considered connected.

    q = [white_squares[0]]
    visited = {white_squares[0]}
    head = 0
    while head < len(q):
        r, c = q[head]
        head += 1
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < SIZE and 0 <= nc < SIZE and \
               grid[nr][nc] == WHITE and (nr, nc) not in visited:
                visited.add((nr, nc))
                q.append((nr, nc))
    
    return len(visited) == len(white_squares)

def is_fully_valid():
    """Performs all final checks on a completed grid."""
    # 1. Check word lengths.
    for r in range(SIZE):
        for c in range(SIZE):
            if grid[r][c] == WHITE:
                # Check for the start of a horizontal word.
                if c == 0 or grid[r][c-1] == BLACK:
                    length = 0
                    for k in range(c, SIZE):
                        if grid[r][k] == WHITE:
                            length += 1
                        else:
                            break
                    if length < 3: return False
                # Check for the start of a vertical word.
                if r == 0 or grid[r-1][c] == BLACK:
                    length = 0
                    for k in range(r, SIZE):
                        if grid[k][c] == WHITE:
                            length += 1
                        else:
                            break
                    if length < 3: return False

    # 2. Check for cheater squares.
    for r in range(SIZE):
        for c in range(SIZE):
            if grid[r][c] == BLACK:
                h_separates = (c > 0 and grid[r][c-1] == WHITE) and \
                              (c < SIZE - 1 and grid[r][c+1] == WHITE)
                v_separates = (r > 0 and grid[r-1][c] == WHITE) and \
                              (r < SIZE - 1 and grid[r+1][c] == WHITE)
                if not (h_separates or v_separates):
                    return False # This is a cheater square.
    
    # 3. Check for connectivity.
    if not check_connectivity():
        return False
        
    return True

def generate_grids(k):
    """Recursively generates and validates grids."""
    global solution_count

    if k == len(cells_to_fill):
        if is_fully_valid():
            solution_count += 1
        return

    r, c = cells_to_fill[k]
    sym_r, sym_c = SIZE - 1 - r, SIZE - 1 - c

    # --- Choice 1: Place WHITE ---
    grid[r][c] = grid[sym_r][sym_c] = WHITE
    generate_grids(k + 1)
    
    # --- Choice 2: Place BLACK ---
    grid[r][c] = grid[sym_r][sym_c] = BLACK
    
    # Pruning: Check for invalid patterns caused by this black square.
    # Check for 2x2 black squares (a common type of cheater).
    # We only need to check the 2x2 block ending at (r,c) due to fill order.
    if r > 0 and c > 0 and \
       grid[r-1][c] == BLACK and grid[r][c-1] == BLACK and grid[r-1][c-1] == BLACK:
       # Backtrack and exit this path.
       grid[r][c] = grid[sym_r][sym_c] = UNDECIDED
       return
    
    # Pruning: Check for newly formed short words.
    # Check horizontal word to the left of (r,c)
    if c > 0 and grid[r][c-1] == WHITE:
        if (c == 1 or grid[r][c-2] == BLACK): # Word is length 1
            grid[r][c] = grid[sym_r][sym_c] = UNDECIDED
            return
        if (c == 2 or grid[r][c-3] == BLACK) and grid[r][c-2] == WHITE: # Word is length 2
            grid[r][c] = grid[sym_r][sym_c] = UNDECIDED
            return

    # Check vertical word above (r,c)
    if r > 0 and grid[r-1][c] == WHITE:
        if (r == 1 or grid[r-2][c] == BLACK): # Word is length 1
            grid[r][c] = grid[sym_r][sym_c] = UNDECIDED
            return
        if (r == 2 or grid[r-3][c] == BLACK) and grid[r-2][c] == WHITE: # Word is length 2
            grid[r][c] = grid[sym_r][sym_c] = UNDECIDED
            return

    generate_grids(k + 1)
    
    # Backtrack to undo the choice for future exploration.
    grid[r][c] = grid[sym_r][sym_c] = UNDECIDED


if __name__ == '__main__':
    # Initial message about the process
    # print("Searching for valid 8x8 crossword grids. This may take some time...")
    
    # Start the recursive search from the first cell.
    generate_grids(0)
    
    # Print the final result.
    # The phrasing "output each number in the final equation" is interpreted as
    # providing the parameters and the final result clearly.
    print(f"For a puzzle of size {SIZE}x{SIZE}, with the given rules, the number of possible grids is:")
    print(solution_count)
<<<2>>>