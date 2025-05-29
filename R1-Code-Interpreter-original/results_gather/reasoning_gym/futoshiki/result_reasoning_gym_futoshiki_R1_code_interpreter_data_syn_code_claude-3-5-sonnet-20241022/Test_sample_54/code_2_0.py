def print_solution(grid):
    result = []
    for i in range(4):
        row = []
        for j in range(4):
            row.append(str(grid[i][j]))
            # Add inequality signs for first row
            if i == 0 and j < 3:
                if j == 0 or j == 2:
                    row.append(" > ")
                else:
                    row.append("   ")
        result.append("   ".join(row))
        if i < 3:
            result.append("             ")
    for line in result:
        print(line)

def is_valid(grid, row, col, num):
    # Check row
    for j in range(4):
        if grid[row][j] == num:
            return False
    
    # Check column
    for i in range(4):
        if grid[i][col] == num:
            return False
            
    # Check given numbers
    if row == 1 and col == 1 and num != 3:
        return False
    if row == 1 and col == 3 and num != 4:
        return False
    if row == 2 and col == 1 and num != 4:
        return False
    if row == 2 and col == 2 and num != 3:
        return False
    if row == 2 and col == 3 and num != 2:
        return False
    
    # Check first row inequalities
    if row == 0:
        # First inequality (position 0 > position 1)
        if col == 0 and grid[0][1] != 0:
            if not (num > grid[0][1]):
                return False
        if col == 1 and grid[0][0] != 0:
            if not (grid[0][0] > num):
                return False
        # Second inequality (position 2 > position 3)
        if col == 2 and grid[0][3] != 0:
            if not (num > grid[0][3]):
                return False
        if col == 3 and grid[0][2] != 0:
            if not (grid[0][2] > num):
                return False
    
    return True

def find_empty(grid):
    for i in range(4):
        for j in range(4):
            if grid[i][j] == 0:
                return (i, j)
    return None

def solve(grid):
    empty = find_empty(grid)
    if not empty:
        # Verify final solution
        for i in range(4):
            # Check rows
            if len(set(grid[i])) != 4:
                return False
            # Check columns
            if len(set(grid[j][i] for j in range(4))) != 4:
                return False
        # Check first row inequalities
        if not (grid[0][0] > grid[0][1] and grid[0][2] > grid[0][3]):
            return False
        return True
    
    row, col = empty
    for num in range(1, 5):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
    return False

# Initialize puzzle
grid = [
    [0, 0, 0, 0],
    [0, 3, 0, 4],
    [0, 4, 3, 2],
    [0, 0, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")