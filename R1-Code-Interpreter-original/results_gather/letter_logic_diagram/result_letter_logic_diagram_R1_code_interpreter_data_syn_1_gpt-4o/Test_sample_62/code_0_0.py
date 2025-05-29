def solve_puzzle(grid):
    from itertools import permutations

    # Check if a letter can be placed at a given position
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
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for letter in 'abcdefg':
                        if can_place(letter, row, col):
                            grid[row][col] = letter
                            if fill_remaining():
                                return True
                            grid[row][col] = ''
                    return False
        return True

    # Start filling the grid
    fill_grid()

    # Return the filled grid
    return grid

# Initial grid with pre-filled cells
initial_grid = [
    ['', 'b', 'a', '', '', '', ''],
    ['b', '', '', 'd', '', '', ''],
    ['', '', 'd', 'f', '', '', 'b'],
    ['g', 'd', '', '', 'e', '', ''],
    ['d', 'f', 'c', 'e', 'b', 'a', ''],
    ['', '', '', '', '', '', 'd'],
    ['', 'e', '', '', '', '', 'f']
]

# Solve the puzzle
solved_grid = solve_puzzle(initial_grid)

# Print the solution
for row in solved_grid:
    print(','.join(row))