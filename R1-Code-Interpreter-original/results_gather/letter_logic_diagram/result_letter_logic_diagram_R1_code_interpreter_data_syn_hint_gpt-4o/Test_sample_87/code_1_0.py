def is_valid(grid, row, col, char):
    # Check if char is not in the current row
    if char in grid[row]:
        return False
    # Check if char is not in the current column
    for r in range(7):
        if grid[r][col] == char:
            return False
    return True

def solve(grid, minor_diagonal_char, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, minor_diagonal_char, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_char, row, col + 1)
    
    # If on the minor diagonal, place the minor diagonal character
    if col == 6 - row:
        if is_valid(grid, row, col, minor_diagonal_char):
            grid[row][col] = minor_diagonal_char
            if solve(grid, minor_diagonal_char, row, col + 1):
                return True
            grid[row][col] = ''
        return False
    
    # Try placing each character
    for char in 'abcdefg':
        if is_valid(grid, row, col, char):
            grid[row][col] = char
            if solve(grid, minor_diagonal_char, row, col + 1):
                return True
            grid[row][col] = ''
    return False

# Initial grid setup
grid = [
    ['e', 'a', 'c', '', '', '', ''],
    ['', '', 'g', 'f', '', 'd', 'e'],
    ['c', 'g', 'f', 'b', 'd', '', ''],
    ['', '', 'b', 'd', 'e', '', 'c'],
    ['', '', '', '', '', '', ''],
    ['', 'd', 'e', '', 'c', 'g', 'f'],
    ['d', '', '', '', 'g', 'f', '']
]

# Try each letter for the minor diagonal
for minor_diagonal_char in 'abcdefg':
    # Make a deep copy of the grid
    grid_copy = [row[:] for row in grid]
    if solve(grid_copy, minor_diagonal_char):
        for row in grid_copy:
            print(','.join(row))
        break