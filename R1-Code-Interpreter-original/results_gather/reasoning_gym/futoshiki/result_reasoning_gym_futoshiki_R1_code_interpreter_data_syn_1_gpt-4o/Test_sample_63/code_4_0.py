def is_valid(grid, row, col, num):
    # Check if num is not in the current row and column
    for i in range(8):
        if grid[row][i] == num or grid[i][col] == num:
            return False

    # Specific constraints from the puzzle
    # Row 1, Column 1 ∧ Row 2, Column 1
    if row == 0 and col == 0 and grid[1][0] != 0 and grid[1][0] <= num:
        return False
    if row == 1 and col == 0 and grid[0][0] != 0 and grid[0][0] >= num:
        return False

    # Row 3, Column 5 ∨ Row 4, Column 5
    if row == 2 and col == 4 and grid[3][4] != 0 and grid[3][4] >= num:
        return False
    if row == 3 and col == 4 and grid[2][4] != 0 and grid[2][4] <= num:
        return False

    # Row 5, Column 6 ∧ Row 6, Column 6
    if row == 4 and col == 5 and grid[5][5] != 0 and grid[5][5] <= num:
        return False
    if row == 5 and col == 5 and grid[4][5] != 0 and grid[4][5] >= num:
        return False

    # Row 7, Column 1 ∧ Row 8, Column 1
    if row == 6 and col == 0 and grid[7][0] != 0 and grid[7][0] <= num:
        return False
    if row == 7 and col == 0 and grid[6][0] != 0 and grid[6][0] >= num:
        return False

    # Row 7, Column 8 ∧ Row 8, Column 8
    if row == 6 and col == 7 and grid[7][7] != 0 and grid[7][7] <= num:
        return False
    if row == 7 and col == 7 and grid[6][7] != 0 and grid[6][7] >= num:
        return False

    # Row 7, Column 6 > Row 7, Column 7
    if row == 6 and col == 5 and grid[6][6] != 0 and grid[6][6] >= num:
        return False
    if row == 6 and col == 6 and grid[6][5] != 0 and grid[6][5] <= num:
        return False

    # Row 8, Column 3 > Row 8, Column 4
    if row == 7 and col == 2 and grid[7][3] != 0 and grid[7][3] >= num:
        return False
    if row == 7 and col == 3 and grid[7][2] != 0 and grid[7][2] <= num:
        return False

    # Row 8, Column 4 > Row 8, Column 5
    if row == 7 and col == 3 and grid[7][4] != 0 and grid[7][4] >= num:
        return False
    if row == 7 and col == 4 and grid[7][3] != 0 and grid[7][3] <= num:
        return False

    # Row 9, Column 1 > Row 9, Column 2
    if row == 8 and col == 0 and grid[8][1] != 0 and grid[8][1] >= num:
        return False
    if row == 8 and col == 1 and grid[8][0] != 0 and grid[8][0] <= num:
        return False

    return True

def solve_futoshiki(grid):
    for row in range(8):
        for col in range(8):
            if grid[row][col] == 0:
                for num in range(1, 9):
                    if is_valid(grid, row, col, num):
                        grid[row][col] = num
                        if solve_futoshiki(grid):
                            return True
                        grid[row][col] = 0
                return False
    return True

# Initial grid with 0 representing empty cells
grid = [
    [0, 0, 3, 0, 7, 0, 1, 8],
    [0, 0, 0, 3, 6, 0, 4, 1],
    [2, 5, 1, 0, 4, 0, 0, 0],
    [1, 2, 0, 0, 0, 4, 0, 3],
    [0, 4, 0, 0, 0, 3, 0, 5],
    [0, 0, 0, 8, 0, 7, 0, 0],
    [0, 8, 5, 0, 0, 6, 0, 0],
    [0, 0, 0, 0, 0, 0, 5, 7]
]

solve_futoshiki(grid)

# Print the solved grid
for row in grid:
    print(row)