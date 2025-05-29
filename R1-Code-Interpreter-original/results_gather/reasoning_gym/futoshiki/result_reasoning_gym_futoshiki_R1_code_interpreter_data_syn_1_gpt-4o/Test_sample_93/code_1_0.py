def is_valid(grid, row, col, num, constraints):
    # Check if num is not in the current row and column
    for i in range(5):
        if grid[row][i] == num or grid[i][col] == num:
            return False
    
    # Check constraints
    for (r1, c1, r2, c2, op) in constraints:
        if (r1 == row and c1 == col) or (r2 == row and c2 == col):
            if op == '>' and grid[r1][c1] != 0 and grid[r2][c2] != 0:
                if not (grid[r1][c1] > grid[r2][c2]):
                    return False
            if op == '<' and grid[r1][c1] != 0 and grid[r2][c2] != 0:
                if not (grid[r1][c1] < grid[r2][c2]):
                    return False
    return True

def solve_futoshiki(grid, constraints, row=0, col=0):
    if row == 5:
        return True
    if col == 5:
        return solve_futoshiki(grid, constraints, row + 1, 0)
    if grid[row][col] != 0:
        return solve_futoshiki(grid, constraints, row, col + 1)
    
    for num in range(1, 6):
        if is_valid(grid, row, col, num, constraints):
            grid[row][col] = num
            if solve_futoshiki(grid, constraints, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

def main():
    constraints = [
        (0, 2, 0, 3, '>'),
        (1, 2, 1, 3, '>'),
        (1, 3, 2, 3, '∨'),
        (2, 1, 3, 1, '∧'),
        (3, 3, 4, 3, '∨')
    ]
    
    grid = [
        [0, 0, 4, 3, 0],
        [1, 0, 0, 2, 0],
        [0, 2, 0, 0, 4],
        [0, 3, 0, 0, 0],
        [0, 0, 0, 0, 0]
    ]
    
    if solve_futoshiki(grid, constraints):
        for row in grid:
            print(row)
    else:
        print("No solution found")

main()