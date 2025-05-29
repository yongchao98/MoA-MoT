def print_solution(grid):
    result = ""
    for i in range(4):
        row = [str(grid[i][j]) if grid[i][j] != 0 else '_' for j in range(4)]
        result += "   ".join(row) + "\n"
        if i < 3:
            if i == 0:
                result += "        ∧   ∨\n"
            elif i == 1:
                result += "∨   ∧        \n"
            elif i == 2:
                result += "             \n"
    print(f"<<<{result}>>>")

def is_valid(grid, row, col, num):
    # Check row
    if num in grid[row]:
        return False
    
    # Check column
    if num in [grid[i][col] for i in range(4)]:
        return False
    
    # Check inequality constraints
    if row == 1 and col == 2:  # cell(2,3) > cell(2,2)
        if grid[1][1] != 0 and num <= grid[1][1]:
            return False
    if row == 1 and col == 1:  # cell(2,2) < cell(2,3)
        if grid[1][2] != 0 and num >= grid[1][2]:
            return False
            
    if row == 2 and col == 2:  # cell(3,3) > cell(3,4)
        if grid[2][3] != 0 and num <= grid[2][3]:
            return False
    if row == 2 and col == 3:  # cell(3,4) < cell(3,3)
        if grid[2][2] != 0 and num >= grid[2][2]:
            return False
            
    if row == 1 and col == 0:  # cell(2,1) < cell(3,1)
        if grid[2][0] != 0 and num >= grid[2][0]:
            return False
    if row == 2 and col == 0:  # cell(3,1) > cell(2,1)
        if grid[1][0] != 0 and num <= grid[1][0]:
            return False
            
    if row == 1 and col == 1:  # cell(2,2) > cell(3,2)
        if grid[2][1] != 0 and num <= grid[2][1]:
            return False
    if row == 2 and col == 1:  # cell(3,2) < cell(2,2)
        if grid[1][1] != 0 and num >= grid[1][1]:
            return False
            
    if row == 0 and col == 2:  # cell(1,3) > cell(2,3)
        if grid[1][2] != 0 and num <= grid[1][2]:
            return False
    if row == 1 and col == 2:  # cell(2,3) < cell(1,3)
        if grid[0][2] != 0 and num >= grid[0][2]:
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
    [0, 2, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 3, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")