from itertools import permutations

# Initial grid with empty cells as None
grid = [
    [None, 'b', 'd', None, 'e', 'g', 'f'],
    [None, None, 'a', 'e', 'g', None, None],
    ['d', 'a', None, None, None, 'c', None],
    ['a', 'e', None, None, None, None, 'b'],
    [None, None, 'f', 'c', 'b', None, 'a'],
    [None, 'f', 'c', None, 'd', None, None],
    ['f', None, 'b', None, None, None, 'g']
]

# Function to check if a letter can be placed at a given position
def can_place(grid, row, col, letter):
    # Check row
    if letter in grid[row]:
        return False
    # Check column
    if letter in [grid[i][col] for i in range(7)]:
        return False
    return True

# Function to solve the grid using backtracking
def solve(grid):
    # Find the minor diagonal letter
    diagonal_positions = [(0, 6), (1, 5), (2, 4), (3, 3), (4, 2), (5, 1), (6, 0)]
    for letter in 'abcdefg':
        if all(can_place(grid, i, j, letter) for i, j in diagonal_positions):
            for i, j in diagonal_positions:
                grid[i][j] = letter
            break

    # Backtracking function to fill the grid
    def backtrack(row, col):
        if row == 7:
            return True
        if col == 7:
            return backtrack(row + 1, 0)
        if grid[row][col] is not None:
            return backtrack(row, col + 1)

        for letter in 'abcdefg':
            if can_place(grid, row, col, letter):
                grid[row][col] = letter
                if backtrack(row, col + 1):
                    return True
                grid[row][col] = None
        return False

    backtrack(0, 0)

solve(grid)

# Print the filled grid
for row in grid:
    print(','.join(letter if letter is not None else '_' for letter in row))