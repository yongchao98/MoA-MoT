def is_valid(grid, row, col, letter):
    # Check if the letter is already in the row
    if letter in grid[row]:
        return False
    # Check if the letter is already in the column
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve_grid(grid):
    # Try to fill the grid using backtracking
    for row in range(7):
        for col in range(7):
            if grid[row][col] is None:
                for letter in 'abcdefg':
                    if is_valid(grid, row, col, letter):
                        grid[row][col] = letter
                        if solve_grid(grid):
                            return True
                        grid[row][col] = None
                return False
    return True

# Initial grid with empty cells as None
grid = [
    [None, None, None, 'e', 'd', 'c', None],
    ['f', None, 'e', None, None, None, None],
    [None, None, 'd', None, None, 'g', None],
    ['e', None, 'c', None, None, 'b', 'f'],
    ['d', None, None, 'b', None, 'a', 'e'],
    ['c', 'g', 'b', None, 'a', 'e', 'd'],
    [None, 'b', 'f', None, None, 'd', None]
]

# Determine the letter for the minor diagonal
for letter in 'abcdefg':
    if all(is_valid(grid, r, 6-r, letter) for r in range(7)):
        for r in range(7):
            grid[r][6-r] = letter
        break

# Solve the grid
solve_grid(grid)

# Print the solved grid
for row in grid:
    print(','.join(letter if letter is not None else ' ' for letter in row))