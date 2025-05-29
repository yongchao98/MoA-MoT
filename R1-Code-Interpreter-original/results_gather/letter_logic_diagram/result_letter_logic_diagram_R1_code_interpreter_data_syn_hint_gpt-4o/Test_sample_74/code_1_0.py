from itertools import permutations

# Initial grid with empty cells represented by None
grid = [
    [None, None, None, None, 'd', 'g', 'c'],
    ['a', None, None, 'd', None, 'c', None],
    ['a', None, 'd', 'g', None, None, None],
    [None, None, 'g', 'c', None, None, 'a'],
    ['d', None, 'c', 'e', None, None, 'b'],
    ['g', 'c', 'e', 'f', 'a', 'b', 'd'],
    ['c', 'e', 'f', 'a', None, 'd', 'g']
]

def is_valid(grid, row, col, letter):
    # Check if the letter can be placed at grid[row][col]
    for i in range(7):
        if grid[row][i] == letter or grid[i][col] == letter:
            return False
    return True

def solve(grid, row=0, col=0):
    if row == 7:
        return True
    if col == 7:
        return solve(grid, row + 1, 0)
    if grid[row][col] is not None:
        return solve(grid, row, col + 1)

    for letter in 'abcdefg':
        if is_valid(grid, row, col, letter):
            grid[row][col] = letter
            if solve(grid, row, col + 1):
                return True
            grid[row][col] = None

    return False

# Try each letter for the minor diagonal
for letter in 'abcdefg':
    # Check if the letter can be placed on the minor diagonal
    if all(grid[i][6-i] in (None, letter) for i in range(7)):
        # Place the letter on the minor diagonal
        for i in range(7):
            grid[i][6-i] = letter
        # Attempt to solve the grid
        if solve(grid):
            break

# Print the filled grid
for row in grid:
    print(','.join(row))