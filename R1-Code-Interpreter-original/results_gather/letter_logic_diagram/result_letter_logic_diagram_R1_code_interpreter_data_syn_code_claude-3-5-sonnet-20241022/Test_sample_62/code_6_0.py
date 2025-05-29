def verify_prefilled(grid):
    prefilled = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]
    for i in range(7):
        for j in range(7):
            if prefilled[i][j] != '' and grid[i][j] != prefilled[i][j]:
                return False
    return True

def get_available_letters(grid, row, col):
    # Get letters already used in row
    row_used = set(grid[row])
    # Get letters already used in column
    col_used = set(grid[i][col] for i in range(7))
    # Available letters are those not used in row or column
    available = set('abcdefg') - row_used - col_used
    return available

def initialize_grid():
    # Start with prefilled positions
    grid = [
        ['', 'b', 'a', '', '', '', ''],
        ['b', '', '', 'd', '', '', ''],
        ['', '', 'd', 'f', '', '', 'b'],
        ['g', 'd', '', '', 'e', '', ''],
        ['d', 'f', 'c', 'e', 'b', 'a', ''],
        ['', '', '', '', '', '', 'd'],
        ['', 'e', '', '', '', '', 'f']
    ]
    return grid

def is_valid_diagonal(grid):
    # Check if all filled positions on minor diagonal are the same
    diagonal = []
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal.append(grid[i][6-i])
    return len(set(diagonal)) <= 1

def solve(grid, pos=0):
    if pos == 49:  # Reached end of grid
        return verify_prefilled(grid)

    row = pos // 7
    col = pos % 7

    # If position is prefilled, move to next position
    if grid[row][col] != '':
        return solve(grid, pos + 1)

    # If on minor diagonal, handle specially
    on_diagonal = (row + col == 6)
    if on_diagonal:
        # Find existing diagonal letter if any
        diag_letter = None
        for i in range(7):
            if grid[i][6-i] != '':
                diag_letter = grid[i][6-i]
                break
        if diag_letter:
            # Must use same letter as rest of diagonal
            if diag_letter in get_available_letters(grid, row, col):
                grid[row][col] = diag_letter
                if solve(grid, pos + 1):
                    return True
            grid[row][col] = ''
            return False

    # Try each available letter
    available = get_available_letters(grid, row, col)
    if on_diagonal and not diag_letter:
        # If on diagonal and no diagonal letter yet, prefer 'c'
        if 'c' in available:
            available = ['c'] + list(available - {'c'})

    for letter in available:
        grid[row][col] = letter
        if is_valid_diagonal(grid) and solve(grid, pos + 1):
            return True
        grid[row][col] = ''

    return False

# Initialize and solve
grid = initialize_grid()
if solve(grid):
    for row in grid:
        print(','.join(row))
else:
    print("No solution found")