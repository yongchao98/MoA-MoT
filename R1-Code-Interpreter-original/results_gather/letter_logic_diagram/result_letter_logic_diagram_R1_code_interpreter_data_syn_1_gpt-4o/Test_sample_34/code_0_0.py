def fill_grid(grid):
    from itertools import permutations

    # Function to check if a letter can be placed at a given position
    def can_place(grid, row, col, letter):
        # Check row
        if letter in grid[row]:
            return False
        # Check column
        for r in range(7):
            if grid[r][col] == letter:
                return False
        return True

    # Try to fill the minor diagonal with a single letter
    for letter in 'abcdefg':
        if all(can_place(grid, r, 6-r, letter) for r in range(7)):
            # Fill the minor diagonal
            for r in range(7):
                grid[r][6-r] = letter
            break

    # Fill the rest of the grid
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                for letter in 'abcdefg':
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        break

    return grid

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
for row in filled_grid:
    print(','.join(row))