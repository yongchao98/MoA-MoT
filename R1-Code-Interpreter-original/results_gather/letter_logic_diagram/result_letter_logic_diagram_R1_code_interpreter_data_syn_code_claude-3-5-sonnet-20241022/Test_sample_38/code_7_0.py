def print_grid(grid):
    for row in grid:
        print(','.join(row))

# Original puzzle - this is our reference
ORIGINAL = [
    ['a','e','','f','g','c',''],
    ['','b','','','','d','a'],
    ['b','f','g','c','d','a','e'],
    ['f','','c','','a','','b'],
    ['g','c','d','a','e','b','f'],
    ['c','','','e','b','','g'],
    ['d','','','','f','g','c']
]

def is_valid_placement(grid, row, col, letter):
    # Check if this conflicts with original puzzle
    if ORIGINAL[row][col] != '' and ORIGINAL[row][col] != letter:
        return False

    # Check row
    for j in range(7):
        if grid[row][j] == letter:
            return False

    # Check column
    for i in range(7):
        if grid[i][col] == letter:
            return False

    # Check minor diagonal (top-right to bottom-left)
    if row + col == 6:  # if on minor diagonal
        for i in range(7):
            j = 6 - i
            if grid[i][j] != '' and grid[i][j] != letter:
                return False

    return True

def solve():
    # Start with empty grid
    grid = [['' for _ in range(7)] for _ in range(7)]
    
    # First, copy all pre-filled values
    for i in range(7):
        for j in range(7):
            if ORIGINAL[i][j] != '':
                grid[i][j] = ORIGINAL[i][j]

    # Find diagonal letter (if any exists)
    diagonal_letter = None
    for i in range(7):
        if grid[i][6-i] != '':
            diagonal_letter = grid[i][6-i]
            break

    def solve_position(pos):
        if pos >= 49:  # 7x7 = 49 cells
            return True

        row = pos // 7
        col = pos % 7

        # If position is pre-filled in original, skip to next position
        if ORIGINAL[row][col] != '':
            return solve_position(pos + 1)

        # If on minor diagonal and we know the diagonal letter
        if row + col == 6 and diagonal_letter:
            if is_valid_placement(grid, row, col, diagonal_letter):
                grid[row][col] = diagonal_letter
                if solve_position(pos + 1):
                    return True
                grid[row][col] = ''
            return False

        # Try each possible letter
        for letter in 'abcdefg':
            if is_valid_placement(grid, row, col, letter):
                grid[row][col] = letter
                if solve_position(pos + 1):
                    return True
                grid[row][col] = ''

        return False

    if solve_position(0):
        print_grid(grid)
        return True
    else:
        print("No solution exists")
        return False

# Solve the puzzle
solve()