def print_solution(grid):
    result = ""
    for i in range(6):
        row = "   ".join(str(grid[i][j]) for j in range(6))
        result += row + "\n"
        if i < 5:
            result += "                     \n"
    print(result)

def is_valid(grid, row, col, num):
    # Check row
    for x in range(6):
        if grid[row][x] == num:
            return False
    
    # Check column
    for x in range(6):
        if grid[x][col] == num:
            return False
    
    # Check initial values match
    if row == 0 and col == 1 and num != 2: return False
    if row == 0 and col == 5 and num != 4: return False
    if row == 1 and col == 0 and num != 6: return False
    if row == 1 and col == 1 and num != 1: return False
    if row == 1 and col == 5 and num != 5: return False
    if row == 2 and col == 1 and num != 6: return False
    if row == 2 and col == 3 and num != 2: return False
    if row == 3 and col == 4 and num != 1: return False
    if row == 4 and col == 3 and num != 3: return False
    if row == 4 and col == 5 and num != 2: return False
    if row == 5 and col == 1 and num != 4: return False
    if row == 5 and col == 2 and num != 6: return False
    if row == 5 and col == 4 and num != 2: return False
    
    # Additional constraints based on puzzle rules
    # Row 1: _ 2 _ _ _ 4 (all numbers before 4 must be less than 4)
    if row == 0 and col < 5 and num >= 4: return False
    
    # Row 2: 6 1 _ _ _ 5 (numbers between 1 and 5 must be in range 2-4)
    if row == 1 and col in [2,3,4] and (num <= 1 or num >= 5): return False
    
    # Row 3: _ 6 _ 2 _ _ (numbers after 6 must be less than 6)
    if row == 2 and col > 1 and num >= 6: return False
    
    # Row 4: _ _ _ _ 1 _ (last number must be greater than 1)
    if row == 3 and col == 5 and num <= 1: return False
    
    # Row 5: _ _ _ 3 _ 2 (number before 2 must be greater than 2)
    if row == 4 and col == 4 and num <= 2: return False
    
    return True

def solve_futoshiki(grid, row=0, col=0):
    if col == 6:
        row += 1
        col = 0
    
    if row == 6:
        return True
    
    if grid[row][col] != 0:
        return solve_futoshiki(grid, row, col + 1)
    
    for num in range(1, 7):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve_futoshiki(grid, row, col + 1):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid (0 represents empty cells)
grid = [
    [0, 2, 0, 0, 0, 4],
    [6, 1, 0, 0, 0, 5],
    [0, 6, 0, 2, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 3, 0, 2],
    [0, 4, 6, 0, 2, 0]
]

if solve_futoshiki(grid):
    print_solution(grid)
else:
    print("No solution exists")