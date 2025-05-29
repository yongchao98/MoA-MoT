def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row or column
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                # Determine if this is a minor diagonal position
                if col == 6 - row:
                    if is_valid(grid, row, col, minor_diagonal_letter):
                        grid[row][col] = minor_diagonal_letter
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[row][col] = ''
                else:
                    for letter in 'abcdefg':
                        if is_valid(grid, row, col, letter):
                            grid[row][col] = letter
                            if solve(grid, minor_diagonal_letter):
                                return True
                            grid[row][col] = ''
                return False
    return True

# Initial grid
grid = [
    ['', 'g', 'a', '', 'd', '', ''],
    ['a', '', '', '', '', 'c', 'f'],
    ['a', '', 'd', 'b', 'c', '', 'g'],
    ['e', 'd', '', 'c', 'f', 'g', 'a'],
    ['', 'b', 'c', '', '', 'a', ''],
    ['', 'c', '', 'g', 'a', 'e', 'd'],
    ['', 'f', '', 'a', '', 'd', 'b']
]

# Determine the letter for the minor diagonal
minor_diagonal_indices = [(i, 6 - i) for i in range(7)]
minor_diagonal_letters = {grid[i][6 - i] for i in range(7) if grid[i][6 - i]}
minor_diagonal_letter = (set('abcdefg') - minor_diagonal_letters).pop()

# Solve the grid
solve(grid, minor_diagonal_letter)

# Print the filled grid
for row in grid:
    print(','.join(row))