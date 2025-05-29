from itertools import permutations

# Initial grid setup
grid = [
    ['e', 'c', 'b', '', 'f', 'd', 'g'],
    ['c', '', '', 'f', '', '', 'e'],
    ['b', 'a', '', '', 'g', 'e', 'c'],
    ['', 'f', '', '', 'e', 'c', ''],
    ['', 'd', '', 'e', '', '', ''],
    ['', '', 'e', 'c', 'b', '', 'f'],
    ['', '', 'c', '', '', '', '']
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

# Backtracking function to fill the grid
def solve(grid):
    # Find the next empty cell
    for r in range(7):
        for c in range(7):
            if grid[r][c] == '':
                # Try each letter
                for letter in 'abcdefg':
                    if can_place(grid, r, c, letter):
                        grid[r][c] = letter
                        if solve(grid):
                            return True
                        grid[r][c] = ''  # Backtrack
                return False
    return True

# Determine the letter for the minor diagonal
for letter in 'abcdefg':
    if all(can_place(grid, r, 6-r, letter) for r in range(7)):
        for r in range(7):
            grid[r][6-r] = letter
        break

# Solve the grid
solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(row))