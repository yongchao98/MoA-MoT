# Define the initial puzzle with constraints
puzzle = [
    [0, 0, 3, 2],
    [0, 0, 0, 0],
    [0, 2, 0, 3],
    [0, 1, 0, 0]
]

# Define the constraints
constraints = {
    'horizontal': [(1, 0, 1)],  # (row, col1, col2) for '<' constraint
    'vertical': [(1, 1, 2)]     # (col, row1, row2) for '∨' constraint
}

# Function to check if a number can be placed in a given position
def is_valid(puzzle, row, col, num):
    # Check row and column uniqueness
    for i in range(4):
        if puzzle[row][i] == num or puzzle[i][col] == num:
            return False
    # Check horizontal constraints
    for (r, c1, c2) in constraints['horizontal']:
        if r == row:
            if c1 == col and puzzle[row][c2] != 0 and num >= puzzle[row][c2]:
                return False
            if c2 == col and puzzle[row][c1] != 0 and num <= puzzle[row][c1]:
                return False
    # Check vertical constraints
    for (c, r1, r2) in constraints['vertical']:
        if c == col:
            if r1 == row and puzzle[r2][col] != 0 and num >= puzzle[r2][col]:
                return False
            if r2 == row and puzzle[r1][col] != 0 and num <= puzzle[r1][col]:
                return False
    return True

# Function to solve the puzzle using backtracking
def solve(puzzle):
    for row in range(4):
        for col in range(4):
            if puzzle[row][col] == 0:
                for num in range(1, 5):
                    if is_valid(puzzle, row, col, num):
                        puzzle[row][col] = num
                        if solve(puzzle):
                            return True
                        puzzle[row][col] = 0
                return False
    return True

# Solve the puzzle
solve(puzzle)

# Print the solved puzzle
for row in puzzle:
    print(row)