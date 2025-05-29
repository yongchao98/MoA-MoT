from itertools import permutations

# Initial grid setup
grid = [
    [' ', 'c', ' ', ' ', ' ', 'd', ' '],
    ['c', ' ', ' ', 'b', 'd', 'f', 'e'],
    ['a', ' ', ' ', ' ', 'f', 'e', ' '],
    [' ', 'b', ' ', ' ', 'e', ' ', ' '],
    [' ', ' ', 'f', ' ', ' ', 'a', ' '],
    ['d', ' ', 'e', ' ', ' ', ' ', ' '],
    [' ', ' ', ' ', 'a', 'g', 'b', ' ']
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
    diagonal_letter = None
    for letter in 'abcdefg':
        if all(grid[i][6-i] in (' ', letter) for i in range(7)):
            diagonal_letter = letter
            break

    # Fill the minor diagonal
    for i in range(7):
        grid[i][6-i] = diagonal_letter

    # Backtracking to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] != ' ':
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if can_place(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = ' '

        return False

    backtrack(0, 0)
    return grid

# Solve the grid
solved_grid = solve(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))