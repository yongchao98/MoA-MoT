def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, row, col + 1)
    
    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = ''
    
    return False

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] in ('', letter) for i in range(7)):
            minor_diagonal_letter = letter
            break
    
    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = minor_diagonal_letter
    
    # Solve the grid
    solve(grid)
    return grid

# Initial grid
grid = [
    ['', 'd', 'c', 'f', '', 'b', 'a'],
    ['', '', 'f', 'e', 'b', 'a', ''],
    ['c', '', 'e', 'b', 'a', '', ''],
    ['', '', 'b', 'a', 'g', 'd', ''],
    ['', '', '', 'g', 'd', 'c', ''],
    ['', '', '', 'd', '', '', 'e'],
    ['', '', '', 'c', 'f', '', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))