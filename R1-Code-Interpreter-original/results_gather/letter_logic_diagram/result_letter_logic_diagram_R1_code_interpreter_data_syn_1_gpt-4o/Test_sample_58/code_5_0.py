def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the given row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid):
    # Try to fill the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                # If it's a diagonal cell, it must be the diagonal letter
                if col == 6 - row:
                    for letter in 'abcdefg':
                        if all(grid[i][6-i] in ('', letter) for i in range(7)):
                            grid[row][col] = letter
                            if solve(grid):
                                return True
                            grid[row][col] = ''
                else:
                    for letter in 'abcdefg':
                        if is_valid(grid, row, col, letter):
                            grid[row][col] = letter
                            if solve(grid):
                                return True
                            grid[row][col] = ''
                return False
    return True

# Initial grid setup
grid = [
    ['d', '', '', '', '', 'f', ''],
    ['', 'b', '', '', 'f', 'a', 'd'],
    ['b', '', '', 'f', '', '', 'd'],
    ['', '', 'f', '', '', '', ''],
    ['c', '', '', '', 'e', '', 'g'],
    ['', 'a', 'd', '', '', '', 'c'],
    ['', '', 'e', '', '', 'c', 'f']
]

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))