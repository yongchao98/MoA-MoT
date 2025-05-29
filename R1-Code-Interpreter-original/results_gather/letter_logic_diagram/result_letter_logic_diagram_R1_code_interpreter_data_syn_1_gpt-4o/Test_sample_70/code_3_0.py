def can_place_on_diagonal(grid, letter):
    for i in range(7):
        if grid[i][6-i] not in ('', letter):
            return False
    return True

def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the current row and column
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

# Initial grid with pre-filled values
grid = [
    ['', '', 'c', 'f', '', '', ''],
    ['e', 'c', 'f', '', '', '', ''],
    ['', 'f', 'g', 'b', 'a', '', 'e'],
    ['f', 'g', 'b', 'a', '', '', 'c'],
    ['g', 'b', 'a', 'd', 'e', '', 'f'],
    ['', '', '', '', 'c', '', ''],
    ['a', '', 'e', '', '', '', '']
]

# Try each letter for the diagonal
for letter in 'abcdefg':
    if can_place_on_diagonal(grid, letter):
        # Place the diagonal letter
        for i in range(7):
            grid[i][6-i] = letter
        if solve(grid):
            break

# Print the filled grid
for row in grid:
    print(','.join(row))