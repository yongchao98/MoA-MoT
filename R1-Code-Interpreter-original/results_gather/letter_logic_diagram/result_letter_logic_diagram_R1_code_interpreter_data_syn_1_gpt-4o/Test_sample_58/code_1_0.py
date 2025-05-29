def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the given row and column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, diagonal_letter):
    # Try to fill the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                # If it's a diagonal cell, it must be the diagonal letter
                if col == 6 - row:
                    grid[row][col] = diagonal_letter
                    if solve(grid, diagonal_letter):
                        return True
                    grid[row][col] = ''
                else:
                    for letter in 'abcdefg':
                        if is_valid(grid, row, col, letter):
                            grid[row][col] = letter
                            if solve(grid, diagonal_letter):
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

# Determine the letter for the minor diagonal
for letter in 'abcdefg':
    if all(grid[i][6-i] in ('', letter) for i in range(7)):
        diagonal_letter = letter
        break

# Solve the grid
solve(grid, diagonal_letter)

# Print the filled grid
for row in grid:
    print(','.join(row))