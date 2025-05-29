def solve_puzzle(grid):
    from collections import Counter

    # Function to check if a letter can be placed at a given position
    def can_place(letter, row, col):
        # Check row
        if letter in grid[row]:
            return False
        # Check column
        for r in range(7):
            if grid[r][col] == letter:
                return False
        return True

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if can_place(letter, row, col):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''

        return False

    # Fill the minor diagonal with a letter that can fit all positions
    def fill_minor_diagonal():
        for letter in 'abcdefg':
            if all(can_place(letter, r, 6-r) for r in range(7)):
                for r in range(7):
                    grid[r][6-r] = letter
                return True
        return False

    # Fill the minor diagonal
    if not fill_minor_diagonal():
        return None

    # Use backtracking to fill the rest of the grid
    if backtrack(0, 0):
        return grid
    else:
        return None

# Initial grid setup
grid = [
    ['', 'c', '', '', '', 'd', ''],
    ['c', '', '', 'b', 'd', 'f', 'e'],
    ['a', '', '', '', 'f', 'e', ''],
    ['', 'b', '', '', 'e', '', ''],
    ['', '', 'f', '', '', 'a', ''],
    ['d', '', 'e', '', '', '', ''],
    ['', '', '', 'a', 'g', 'b', '']
]

# Solve the puzzle
solved_grid = solve_puzzle(grid)

# Print the solved grid
if solved_grid:
    for row in solved_grid:
        print(','.join(row))
else:
    print("No solution found")