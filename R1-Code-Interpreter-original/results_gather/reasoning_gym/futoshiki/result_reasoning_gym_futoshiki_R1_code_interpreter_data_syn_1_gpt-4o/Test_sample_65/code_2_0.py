# Define the constraints
horizontal_constraints = [
    (0, 0, '<'), (0, 1, '>'), (0, 2, '<'),
    (1, 0, '<'), (1, 1, '<'), (1, 2, '>'),
    (2, 2, '<'), (2, 3, '<'),
    (3, 1, '>'),
    (4, 0, '>'),
]

vertical_constraints = [
    (0, 0, '∧'),
    (1, 2, '∨'), (1, 4, '∧'),
    (3, 3, '∨'),
]

def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Check horizontal constraints
    if col < 4:
        if (row, col, '<') in horizontal_constraints and not (num < grid[row][col + 1]):
            return False
        if (row, col, '>') in horizontal_constraints and not (num > grid[row][col + 1]):
            return False
    if col > 0:
        if (row, col - 1, '<') in horizontal_constraints and not (grid[row][col - 1] < num):
            return False
        if (row, col - 1, '>') in horizontal_constraints and not (grid[row][col - 1] > num):
            return False

    # Check vertical constraints
    if row < 4:
        if (row, col, '∧') in vertical_constraints and not (num > grid[row + 1][col]):
            return False
        if (row, col, '∨') in vertical_constraints and not (num < grid[row + 1][col]):
            return False
    if row > 0:
        if (row - 1, col, '∧') in vertical_constraints and not (grid[row - 1][col] > num):
            return False
        if (row - 1, col, '∨') in vertical_constraints and not (grid[row - 1][col] < num):
            return False

    return True

def solve_futoshiki(grid, row=0, col=0):
    if row == 5:
        return True
    if col == 5:
        return solve_futoshiki(grid, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)

    for num in range(1, 6):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0

    return False

# Initialize the grid with zeros
grid = [[0] * 5 for _ in range(5)]

# Solve the puzzle
if solve_futoshiki(grid):
    for row in range(5):
        line = ' '.join(str(grid[row][col]) for col in range(5))
        if row < 4:
            line += ' ' + ' '.join('<' if (row, col, '<') in horizontal_constraints else '>' if (row, col, '>') in horizontal_constraints else ' ' for col in range(4))
        print(line)
        if row < 4:
            line = ' '.join('∧' if (row, col, '∧') in vertical_constraints else '∨' if (row, col, '∨') in vertical_constraints else ' ' for col in range(5))
            print(line)
else:
    print("No solution found")