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
    
    for char in 'abcdefg':
        if (row + col == 6 and char != minor_diagonal_char) or not is_valid(grid, row, col, char):
            continue
        grid[row][col] = char
        if solve(grid, minor_diagonal_char, row, col + 1):
            return True
        grid[row][col] = ''
    return False

def fill_grid(grid):
    # Try each letter for the minor diagonal
    for char in 'abcdefg':
        temp_grid = [row[:] for row in grid]
        for i in range(7):
            temp_grid[i][6-i] = char
        if solve(temp_grid, char):
            return temp_grid
    return None

# Initial grid setup
grid = [
    ['g', 'b', 'e', '', 'd', 'f', ''],
    ['b', '', 'c', '', '', '', 'a'],
    ['e', 'c', '', 'f', '', '', 'b'],
    ['', 'd', 'f', '', '', '', ''],
    ['', '', '', '', 'b', 'e', 'c'],
    ['', '', '', '', '', '', 'd'],
    ['', '', 'b', '', '', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")