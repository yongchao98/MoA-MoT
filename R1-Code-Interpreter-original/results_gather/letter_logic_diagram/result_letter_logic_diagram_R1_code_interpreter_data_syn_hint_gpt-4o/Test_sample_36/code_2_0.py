def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, minor_diagonal_letter, row + 1, 0)
    if grid[row][col] != '':
        return solve(grid, minor_diagonal_letter, row, col + 1)

    for letter in 'abcdefg':
        if (row + col == 6 and letter != minor_diagonal_letter) or not is_valid(grid, row, col, letter):
            continue
        grid[row][col] = letter
        if solve(grid, minor_diagonal_letter, row, col + 1):
            return True
        grid[row][col] = ''

    return False

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_indices = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    possible_letters = set('abcdefg')

    for i, j in minor_diagonal_indices:
        if grid[i][j] in possible_letters:
            possible_letters.remove(grid[i][j])

    for minor_diagonal_letter in possible_letters:
        # Try each possible letter for the minor diagonal
        temp_grid = [row[:] for row in grid]
        for i, j in minor_diagonal_indices:
            temp_grid[i][j] = minor_diagonal_letter

        # Solve the grid
        if solve(temp_grid, minor_diagonal_letter):
            return temp_grid

    return None

# Initial grid setup
grid = [
    ['', '', '', 'e', 'd', 'c', ''],
    ['f', '', 'e', '', '', '', ''],
    ['', '', 'd', '', '', 'g', ''],
    ['e', '', 'c', '', '', 'b', 'f'],
    ['d', '', '', 'b', '', 'a', 'e'],
    ['c', 'g', 'b', '', 'a', 'e', 'd'],
    ['', 'b', 'f', '', '', 'd', '']
]

filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")