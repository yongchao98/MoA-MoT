from itertools import permutations

# Initial grid with empty cells represented as ''
grid = [
    ['', 'f', 'b', 'e', 'g', 'd', 'a'],
    ['', 'b', 'e', 'g', '', 'a', 'c'],
    ['', '', '', 'd', 'a', '', ''],
    ['e', '', '', 'a', 'c', 'f', ''],
    ['g', '', 'a', 'c', 'f', 'b', ''],
    ['d', 'a', 'c', 'f', 'b', 'e', 'g'],
    ['a', 'c', '', '', 'e', 'g', 'd']
]

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

# Function to solve the grid
def solve(grid):
    # Determine the letter for the minor diagonal
    minor_diagonal_letter = None
    for r, c in zip(range(7), reversed(range(7))):
        if grid[r][c] != '':
            minor_diagonal_letter = grid[r][c]
            break
    if minor_diagonal_letter is None:
        minor_diagonal_letter = 'a'  # Default to 'a' if none is pre-filled

    # Fill the minor diagonal
    for r, c in zip(range(7), reversed(range(7))):
        grid[r][c] = minor_diagonal_letter

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != '':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if can_place(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ''
        return False

    backtrack(0, 0)
    return grid

# Solve the grid
solved_grid = solve(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))