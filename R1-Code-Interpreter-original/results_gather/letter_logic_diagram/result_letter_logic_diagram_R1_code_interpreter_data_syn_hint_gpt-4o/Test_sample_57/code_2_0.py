def is_valid(grid, r, c, letter):
    # Check if the letter can be placed at grid[r][c]
    for i in range(7):
        if grid[r][i] == letter or grid[i][c] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    # Try to fill the grid using backtracking
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                for letter in 'abcdefg':
                    if is_valid(grid, r, c, letter):
                        grid[r][c] = letter
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[r][c] = ''
                return False
    return True

def fill_grid(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    for letter in 'abcdefg':
        if all(is_valid(grid, r, c, letter) for r, c in minor_diagonal_positions):
            for r, c in minor_diagonal_positions:
                grid[r][c] = letter
            if solve(grid, letter):
                return grid
            for r, c in minor_diagonal_positions:
                grid[r][c] = ''
    return None

# Initial grid
grid = [
    ['', '', 'd', 'e', 'a', '', ''],
    ['', '', '', '', '', '', ''],
    ['', '', 'a', '', '', 'f', 'g'],
    ['', '', '', 'b', '', 'g', 'd'],
    ['', 'c', '', '', '', 'd', 'e'],
    ['c', '', 'f', '', '', '', 'a'],
    ['', 'f', 'g', '', 'e', '', 'c']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("No solution found")