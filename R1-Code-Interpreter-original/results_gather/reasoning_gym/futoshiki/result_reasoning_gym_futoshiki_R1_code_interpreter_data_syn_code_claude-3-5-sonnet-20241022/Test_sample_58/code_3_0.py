def print_solution(grid):
    for i in range(5):
        row = []
        for j in range(5):
            row.append(str(grid[i][j]))
        print("   ".join(row))

def check_constraints(grid, row, col, num):
    # Check row
    for j in range(5):
        if j != col and grid[row][j] == num:
            return False
    
    # Check column
    for i in range(5):
        if i != row and grid[i][col] == num:
            return False
    
    # Check vertical constraints
    if col == 1:  # Column 2
        if row == 1:  # Row 2
            # Must be greater than the value below it
            if grid[2][1] != 0 and num <= grid[2][1]:
                return False
        if row == 2:  # Row 3
            # Must be less than the value above it
            if grid[1][1] != 0 and num >= grid[1][1]:
                return False
    
    if col == 4:  # Column 5
        if row == 0:  # Row 1
            # Must be greater than the value below it
            if grid[1][4] != 0 and num <= grid[1][4]:
                return False
        if row == 1:  # Row 2
            # Must be less than the value above it
            if grid[0][4] != 0 and num >= grid[0][4]:
                return False
    
    return True

def solve(grid):
    # Find empty location
    empty = None
    min_possibilities = 6
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                # Count valid possibilities for this cell
                count = 0
                possible_nums = []
                for num in range(1, 6):
                    if check_constraints(grid, i, j, num):
                        count += 1
                        possible_nums.append(num)
                if count < min_possibilities:
                    min_possibilities = count
                    empty = (i, j)
                    best_nums = possible_nums
                if count == 0:  # If any cell has no valid numbers, fail fast
                    return False
    
    if empty is None:
        return True
    
    row, col = empty
    
    # Try the possible numbers for the cell with fewest options
    for num in best_nums:
        if check_constraints(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    
    return False

# Initial grid
grid = [
    [0, 0, 0, 3, 4],
    [5, 0, 0, 4, 0],
    [0, 0, 3, 0, 0],
    [0, 5, 0, 0, 1],
    [0, 0, 0, 0, 2]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")