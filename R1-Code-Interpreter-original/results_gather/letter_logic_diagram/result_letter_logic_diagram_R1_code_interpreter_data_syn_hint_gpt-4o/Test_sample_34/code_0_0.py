def solve_puzzle(grid):
    from itertools import permutations

    # Check if a letter can be placed at a specific position
    def can_place(letter, row, col):
        return letter not in grid[row] and all(grid[r][col] != letter for r in range(7))

    # Try to fill the grid
    def fill_grid():
        # Try each letter for the minor diagonal
        for letter in 'abcdefg':
            # Check if the letter can be placed on the minor diagonal
            if all(can_place(letter, r, 6-r) for r in range(7)):
                # Place the letter on the minor diagonal
                for r in range(7):
                    grid[r][6-r] = letter

                # Try to fill the rest of the grid
                if fill_remaining():
                    return True

                # Reset the minor diagonal
                for r in range(7):
                    grid[r][6-r] = ''

        return False

    # Fill the remaining cells
    def fill_remaining():
        for r in range(7):
            for c in range(7):
                if grid[r][c] == '':
                    for letter in 'abcdefg':
                        if can_place(letter, r, c):
                            grid[r][c] = letter
                            if fill_remaining():
                                return True
                            grid[r][c] = ''
                    return False
        return True

    # Start filling the grid
    fill_grid()
    return grid

# Initial grid setup
grid = [
    ['f', '', '', 'c', 'a', '', 'g'],
    ['e', 'b', 'c', 'a', 'd', 'g', ''],
    ['', 'c', '', '', 'g', 'f', 'e'],
    ['c', 'a', '', '', '', 'e', 'b'],
    ['a', 'd', 'g', '', 'e', '', 'c'],
    ['d', 'g', '', 'e', '', 'c', ''],
    ['', '', '', 'b', '', 'a', '']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))