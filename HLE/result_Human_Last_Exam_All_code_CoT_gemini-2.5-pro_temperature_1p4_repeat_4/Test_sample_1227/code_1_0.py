import sys

def solve_crossword_grids():
    """
    This function calculates the number of valid 8x8 crossword grids based on a set of rules.
    It encapsulates the entire logic and variables to avoid global scope issues.
    """
    
    # The grid is a 2D list where 0 represents a white square and 1 represents a black square.
    grid = [[0] * 8 for _ in range(8)]
    
    # `cells_to_set` will hold the 18 representative cells from the inner 6x6 grid.
    # We only need to decide the color for these cells; the symmetry rule determines the rest.
    cells_to_set = []
    for r in range(1, 4):      # Rows 1, 2, and 3 are in the top half of the grid
        for c in range(1, 7):  # Columns 1 through 6 define the inner grid area
            cells_to_set.append((r,c))

    def check_word_length(g):
        """Checks if all words are at least 3 letters long."""
        # Check rows
        for r in range(1, 7):
            c = 1
            while c < 7:
                if g[r][c] == 0 and g[r][c-1] == 1: # Start of a word
                    length = 0
                    k = c
                    while k < 7 and g[r][k] == 0:
                        length += 1
                        k += 1
                    if g[r][k] == 1 and length < 3:
                        return False
                    c = k
                else:
                    c += 1
        # Check columns
        for c in range(1, 7):
            r = 1
            while r < 7:
                if g[r][c] == 0 and g[r-1][c] == 1: # Start of a word
                    length = 0
                    k = r
                    while k < 7 and g[k][c] == 0:
                        length += 1
                        k += 1
                    if g[k][c] == 1 and length < 3:
                        return False
                    r = k
                else:
                    r += 1
        return True

    def check_connectivity(g):
        """Checks if all white squares are connected using Breadth-First Search (BFS)."""
        first_white, total_white = None, 0
        for r in range(1, 7):
            for c in range(1, 7):
                if g[r][c] == 0:
                    if not first_white:
                        first_white = (r, c)
                    total_white += 1
        
        if total_white == 0:
            return False # A valid crossword must have at least one word.

        q, visited = [first_white], {first_white}
        head = 0
        while head < len(q):
            r, c = q[head]; head += 1
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if g[nr][nc] == 0 and (nr, nc) not in visited:
                    visited.add((nr, nc)); q.append((nr, nc))

        return len(visited) == total_white

    def check_cheaters(g):
        """Checks that there are no 'cheater' (unnecessary) black squares."""
        for r in range(1, 7):
            for c in range(1, 7):
                if g[r][c] == 1: # For each black square
                    # Check if it separates words horizontally or vertically
                    h_sep = g[r][c-1] == 0 and g[r][c+1] == 0
                    v_sep = g[r-1][c] == 0 and g[r+1][c] == 0
                    if not h_sep and not v_sep:
                        return False # Found a cheater
        return True

    # Use a mutable list to hold the count, avoiding complex global/nonlocal state.
    count_ref = [0]
    
    def recursive_search(k):
        """The backtracking function to explore all possible grid patterns."""
        if k == len(cells_to_set): # Base case: grid is fully colored
            if check_word_length(grid) and check_connectivity(grid) and check_cheaters(grid):
                count_ref[0] += 1
            return

        r, c = cells_to_set[k]
        r_sym, c_sym = 7 - r, 7 - c
        
        # Trial 1: Set the symmetric pair to BLACK and recurse
        grid[r][c], grid[r_sym][c_sym] = 1, 1
        recursive_search(k + 1)

        # Trial 2: Set the symmetric pair to WHITE and recurse
        grid[r][c], grid[r_sym][c_sym] = 0, 0
        recursive_search(k + 1)

    # --- Main execution ---
    # 1. Initialize the grid with a black perimeter and a white interior.
    for i in range(8):
        grid[0][i] = grid[7][i] = 1
        grid[i][0] = grid[i][7] = 1
    for r in range(1, 7):
        for c in range(1, 7):
            grid[r][c] = 0
            
    # 2. Start the search from the first pair of cells.
    recursive_search(0)
    final_count = count_ref[0]
    
    print("The final calculation is based on the following numbers:")
    print("Side length of inner grid = 8 - 2 = 6")
    print("Number of independent symmetric pairs = (6 * 6) / 2 = 18")
    print("Total symmetric patterns to test = 2^18 = 262144")
    print(f"Number of patterns passing all checks = {final_count}")

# Execute the main function to find and print the solution.
solve_crossword_grids()
<<<76>>>