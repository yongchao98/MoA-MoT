from itertools import permutations

# Initial grid setup
grid = [
    ['g', '', '', 'c', '', '', 'd'],
    ['', 'b', '', 'e', '', 'd', ''],
    ['b', '', 'e', 'f', 'd', '', 'a'],
    ['c', 'e', 'f', '', '', 'a', 'b'],
    ['e', 'f', 'd', 'g', 'a', 'b', ''],
    ['f', '', '', 'a', 'b', '', ''],
    ['', 'g', 'a', 'b', '', '', 'f']
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
    # Try to fill the minor diagonal with the same letter
    for letter in 'abcdefg':
        # Check if the letter can be placed on the minor diagonal
        if all(can_place(grid, i, 6-i, letter) for i in range(7)):
            # Place the letter on the minor diagonal
            for i in range(7):
                grid[i][6-i] = letter
            break

    # Fill the rest of the grid
    for row in range(7):
        for col in range(7):
            if grid[row][col] == '':
                for letter in 'abcdefg':
                    if can_place(grid, row, col, letter):
                        grid[row][col] = letter
                        break

    return grid

# Solve the grid
solved_grid = solve(grid)

# Print the solved grid
for row in solved_grid:
    print(','.join(row))