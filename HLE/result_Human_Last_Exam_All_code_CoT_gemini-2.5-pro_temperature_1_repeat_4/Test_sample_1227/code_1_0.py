import collections

def solve_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules.

    The rules are:
    1. 180-degree rotational symmetry.
    2. A minimum word length of 3 letters.
    3. All white squares are fully interconnected.
    4. A standard (assumed) black border around the grid.
    5. No "cheater" squares (interpreted as being enforced by the min word length rule).
    """

    # The independent cells are the top half of the inner 6x6 grid.
    # Symmetry is (r, c) <-> (7-r, 7-c).
    # We choose the cells in the top-left portion of the grid to iterate through.
    # Rows 1 & 2: (1,1)...(1,6) and (2,1)...(2,6) -> 12 cells
    # Row 3: (3,1)...(3,6) -> 6 cells
    # Total independent cells to decide = 12 + 6 = 18.
    independent_cells = []
    for r in range(1, 4):
        for c in range(1, 7):
            independent_cells.append((r, c))

    valid_grid_count = 0
    num_possibilities = 2**len(independent_cells)

    for i in range(num_possibilities):
        # 0 represents a white square, 1 represents a black square.
        # Start with a grid that has black borders and is white inside.
        grid = [[0] * 8 for _ in range(8)]
        for r_idx in range(8):
            grid[r_idx][0] = 1
            grid[r_idx][7] = 1
        for c_idx in range(8):
            grid[0][c_idx] = 1
            grid[7][c_idx] = 1

        # Use the bits of 'i' to determine the color of each independent cell.
        temp_i = i
        for r, c in independent_cells:
            if (temp_i & 1) == 1:
                # Set this cell and its symmetric counterpart to black.
                grid[r][c] = 1
                grid[7 - r][7 - c] = 1
            temp_i >>= 1
        
        # Now, validate the fully constructed grid.
        is_len_valid, word_count = check_word_lengths(grid)
        
        # A valid grid must have all words of length >= 3, at least one word,
        # and all white squares must be connected.
        if is_len_valid and word_count > 0:
            if check_connectivity(grid):
                valid_grid_count += 1
    
    print(f"Total possibilities checked for the inner grid: 2^18 = {num_possibilities}")
    print(f"Final count of valid 8x8 crossword grids: {valid_grid_count}")


def check_word_lengths(grid):
    """
    Checks if all words in the grid are of length 3 or more.
    Returns a tuple: (is_valid, number_of_words).
    """
    word_count = 0
    size = 8
    # Check rows for words
    for r in range(size):
        for c in range(size):
            # Find the start of a horizontal word (white square preceded by black or border)
            if grid[r][c] == 0 and (c == 0 or grid[r][c-1] == 1):
                length = 0
                for k in range(c, size):
                    if grid[r][k] == 0:
                        length += 1
                    else:
                        break  # Word ends at a black square
                if length > 0:
                    if length < 3:
                        return (False, 0)
                    word_count += 1

    # Check columns for words
    for c in range(size):
        for r in range(size):
            # Find the start of a vertical word
            if grid[r][c] == 0 and (r == 0 or grid[r-1][c] == 1):
                length = 0
                for k in range(r, size):
                    if grid[k][c] == 0:
                        length += 1
                    else:
                        break  # Word ends at a black square
                if length > 0:
                    if length < 3:
                        return (False, 0)
                    word_count += 1
    
    return (True, word_count)

def check_connectivity(grid):
    """
    Checks if all white squares in the grid form a single connected component.
    Returns True if connected, False otherwise.
    """
    white_squares = []
    for r in range(8):
        for c in range(8):
            if grid[r][c] == 0:
                white_squares.append((r, c))

    if not white_squares:
        return True  # An all-black grid is trivially connected.

    total_white_count = len(white_squares)
    
    # BFS to find all connected white squares from the first one.
    q = collections.deque([white_squares[0]])
    visited = {white_squares[0]}
    
    while q:
        r, c = q.popleft()
        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < 8 and 0 <= nc < 8 and grid[nr][nc] == 0 and (nr, nc) not in visited:
                visited.add((nr, nc))
                q.append((nr, nc))
    
    return len(visited) == total_white_count

if __name__ == '__main__':
    solve_crossword_grids()