def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if (row == 6 - col and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
                        continue
                    grid[row][col] = letter
                    if solve(grid, minor_diagonal_letter):
                        return True
                    grid[row][col] = ''
                return False
    return True

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for i in range(7):
        if grid[i][6 - i] != '':
            if minor_diagonal_letter is None:
                minor_diagonal_letter = grid[i][6 - i]
            elif minor_diagonal_letter != grid[i][6 - i]:
                raise ValueError("Conflicting letters on the minor diagonal")

    # If no letter is pre-filled on the minor diagonal, choose one
    if minor_diagonal_letter is None:
        minor_diagonal_letter = 'a'  # Arbitrarily choose 'a'

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6 - i] = minor_diagonal_letter

    # Solve the grid
    if not solve(grid, minor_diagonal_letter):
        raise ValueError("No solution exists")

    return grid

# Initial grid setup
grid = [
    ['', 'g', 'f', '', 'e', '', ''],
    ['f', 'a', '', 'c', '', '', 'b'],
    ['', '', '', 'c', 'd', '', ''],
    ['', '', 'c', 'd', 'b', '', ''],
    ['e', 'c', '', '', '', '', 'a'],
    ['', 'd', '', 'g', 'f', 'a', 'e'],
    ['', '', 'g', 'f', '', 'e', 'c']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
for row in filled_grid:
    print(','.join(row))