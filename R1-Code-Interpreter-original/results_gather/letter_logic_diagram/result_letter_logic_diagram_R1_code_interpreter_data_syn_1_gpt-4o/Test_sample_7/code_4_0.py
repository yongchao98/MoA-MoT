def solve_puzzle(grid):
    # Function to check if a letter can be placed at a given position
    def can_place(letter, row, col):
        return letter not in grid[row] and all(grid[r][col] != letter for r in range(7))

    # Function to fill the minor diagonal with a given letter
    def fill_minor_diagonal(grid, letter):
        for i in range(7):
            grid[i][6-i] = letter

    # Backtracking function to fill the grid
    def fill_grid(grid):
        for row in range(7):
            for col in range(7):
                if grid[row][col] == '':
                    for letter in 'abcdefg':
                        if can_place(letter, row, col):
                            grid[row][col] = letter
                            if fill_grid(grid):
                                return True
                            grid[row][col] = ''
                    return False
        return True

    # Try each letter for the minor diagonal
    for letter in 'abcdefg':
        # Make a copy of the grid to try this letter
        temp_grid = [row[:] for row in grid]
        fill_minor_diagonal(temp_grid, letter)

        # Check if this letter can be placed on the minor diagonal
        if all(can_place(letter, i, 6-i) for i in range(7)):
            # Try to fill the rest of the grid
            if fill_grid(temp_grid):
                return temp_grid

    return None

# Initial grid
grid = [
    ['', 'g', 'f', '', 'e', '', ''],
    ['', 'f', 'a', '', 'c', '', 'b'],
    ['', '', '', 'c', 'd', '', ''],
    ['', '', 'c', 'd', 'b', '', ''],
    ['e', 'c', '', '', '', '', 'a'],
    ['', 'd', '', 'g', 'f', 'a', 'e'],
    ['', '', 'g', 'f', '', 'e', 'c']
]

# Solve the puzzle
solution = solve_puzzle(grid)

# Print the solution
if solution:
    for row in solution:
        print(','.join(row))
else:
    print("No solution found")