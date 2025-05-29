def print_solution(grid):
    result = ""
    for i in range(4):
        row = [str(grid[i][j]) if grid[i][j] != 0 else "_" for j in range(4)]
        result += "   ".join(row) + "\n"
        if i < 3:
            result += "             \n"
    return result

def is_valid(grid, row, col, num):
    # Check if position should be fixed
    if (row == 0 and col == 0 and num != 4) or \
       (row == 1 and col == 0 and num != 3) or \
       (row == 1 and col == 3 and num != 2) or \
       (row == 3 and col == 2 and num != 3):
        return False
    
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(4)]:
        return False
    
    # Check vertical constraints
    if row == 2 and col == 1:  # Third row, second column (>)
        if grid[3][1] != 0 and num <= grid[3][1]:
            return False
    if row == 3 and col == 1:  # Fourth row, second column
        if grid[2][1] != 0 and grid[2][1] <= num:
            return False
            
    if row == 2 and col == 2:  # Third row, third column (<)
        if grid[3][2] != 0 and num >= grid[3][2]:
            return False
    if row == 3 and col == 2:  # Fourth row, third column
        if grid[2][2] != 0 and grid[2][2] >= num:
            return False
    
    return True

def solve(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                for num in range(1, 5):
                    if is_valid(grid, i, j, num):
                        grid[i][j] = num
                        if solve(grid):
                            return True
                        grid[i][j] = 0
                return False
    return True

# Initialize grid with known values
grid = [
    [4, 0, 0, 0],
    [3, 0, 0, 2],
    [0, 0, 0, 0],
    [0, 0, 3, 0]
]

if solve(grid):
    solution = print_solution(grid)
    print(solution)
else:
    print("No solution exists")