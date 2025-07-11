import collections

def solve_crossword_grids():
    """
    Calculates the number of valid 8x8 crossword grids based on a set of rules.

    The rules are:
    1. 8x8 grid with black borders.
    2. 180-degree rotational symmetry.
    3. Minimum word length of 3 (no 1- or 2-letter words).
    4. All white squares must be interconnected.
    5. No "cheater" squares (inner black squares must touch words on both axes).
    """
    SIZE = 8
    
    # grid: 0 for black, 1 for white. Initialize with black borders.
    grid = [[0] * SIZE for _ in range(SIZE)]
    
    # We only need to decide the colors for the inner 6x6 grid.
    # Due to symmetry, we only iterate through the top half of these cells.
    # For an 8x8 grid, these are the cells in rows 1, 2, 3 and cols 1-6.
    independent_cells = []
    for r in range(1, SIZE // 2):
        for c in range(1, SIZE - 1):
            independent_cells.append((r, c))

    # This will be our final count
    count = 0

    def check_word_length():
        # Check for horizontal words with length < 3
        for r in range(SIZE):
            length = 0
            for c in range(SIZE + 1):
                is_white = c < SIZE and grid[r][c] == 1
                if is_white:
                    length += 1
                else:
                    if 0 < length < 3:
                        return False
                    length = 0
        
        # Check for vertical words with length < 3
        for c in range(SIZE):
            length = 0
            for r in range(SIZE + 1):
                is_white = r < SIZE and grid[r][c] == 1
                if is_white:
                    length += 1
                else:
                    if 0 < length < 3:
                        return False
                    length = 0
        return True

    def check_connectivity():
        white_squares = []
        for r in range(SIZE):
            for c in range(SIZE):
                if grid[r][c] == 1:
                    white_squares.append((r, c))

        if not white_squares:
            return False  # A valid puzzle must have words.

        # BFS to check if all white squares are connected
        q = collections.deque([white_squares[0]])
        visited = {white_squares[0]}
        
        while q:
            r, c = q.popleft()
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                nr, nc = r + dr, c + dc
                if (0 <= nr < SIZE and 0 <= nc < SIZE) and \
                   grid[nr][nc] == 1 and (nr, nc) not in visited:
                    visited.add((nr, nc))
                    q.append((nr, nc))
        
        return len(visited) == len(white_squares)

    def check_no_cheaters():
        # An inner black square must touch white squares on both horizontal and vertical axes.
        for r in range(1, SIZE - 1):
            for c in range(1, SIZE - 1):
                if grid[r][c] == 0:
                    has_horizontal_white = (grid[r][c-1] == 1 or grid[r][c+1] == 1)
                    has_vertical_white = (grid[r-1][c] == 1 or grid[r+1][c] == 1)
                    if not (has_horizontal_white and has_vertical_white):
                        return False
        return True

    def is_valid_grid():
        if not check_word_length(): return False
        if not check_no_cheaters(): return False
        if not check_connectivity(): return False
        return True

    def find_patterns(cell_index=0):
        nonlocal count
        if cell_index == len(independent_cells):
            if is_valid_grid():
                count += 1
            return

        r, c = independent_cells[cell_index]
        r_sym, c_sym = SIZE - 1 - r, SIZE - 1 - c

        # Branch 1: Set the pair of squares to black
        grid[r][c] = 0
        grid[r_sym][c_sym] = 0
        find_patterns(cell_index + 1)

        # Branch 2: Set the pair of squares to white
        grid[r][c] = 1
        grid[r_sym][c_sym] = 1
        find_patterns(cell_index + 1)

    find_patterns()
    print(f"The number of possible 8x8 crossword grids is: {count}")

if __name__ == '__main__':
    solve_crossword_grids()