def is_valid(grid, row, col, letter):
    # Check if the letter can be placed in the given row and column
    if letter in grid[row]:
        return False
    for r in range(7):
        if grid[r][col] == letter:
            return False
    return True

def solve(grid, minor_diagonal_letter):
    # Find the next empty cell
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                # Try placing each letter
                for letter in 'abcdefg':
                    if is_valid(grid, r, c, letter):
                        grid[r][c] = letter
                        # Check if the minor diagonal constraint is satisfied
                        if r + c == 6 and letter != minor_diagonal_letter:
                            grid[r][c] = ''
                            continue
                        if solve(grid, minor_diagonal_letter):
                            return True
                        grid[r][c] = ''
                return False
    return True

def fill_grid(grid):
    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        # Create a copy of the grid to try this letter
        grid_copy = [row[:] for row in grid]
        # Fill the minor diagonal with this letter
        for r in range(7):
            grid_copy[r][6-r] = letter
        # Try to solve the grid
        if solve(grid_copy, letter):
            return grid_copy
    return None

# Initial grid
grid = [
    ['f', '', '', 'c', 'a', '', 'g'],
    ['e', 'b', 'c', 'a', 'd', 'g', ''],
    ['', 'c', '', '', 'g', 'f', 'e'],
    ['c', 'a', '', '', '', 'e', 'b'],
    ['a', 'd', 'g', '', 'e', '', 'c'],
    ['d', 'g', '', 'e', '', 'c', ''],
    ['', '', '', 'b', '', 'a', '']
]

# Fill the grid
filled_grid = fill_grid(grid)

# Print the filled grid
if filled_grid:
    for row in filled_grid:
        print(','.join(row))
else:
    print("Failed to fill the grid correctly.")