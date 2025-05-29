def print_puzzle(grid):
    result = []
    # Horizontal constraints (row, col pairs where < appears)
    h_cons = [(1,6), (3,6), (4,7), (7,4)]
    # Vertical constraints (row, col, type pairs where v=down, ^=up)
    v_cons = [(0,1,'v'), (1,1,'^'), (1,7,'^'), (2,7,'v'), (4,4,'^'), (5,3,'v')]
    
    for i in range(8):
        row = ""
        for j in range(8):
            row += str(grid[i][j])
            if j < 7:
                if (i,j) in h_cons:
                    row += " < "
                else:
                    row += "   "
        result.append(row)
        if i < 7:
            v_row = ["    "] * 8
            for r,c,t in v_cons:
                if r == i:
                    v_row[c] = " ∨  " if t == 'v' else " ∧  "
            result.append("".join(v_row).rstrip())
    
    print("<<<")
    for line in result:
        print(line)
    print(">>>")

def is_valid(grid, row, col, num):
    # Check row
    for x in range(8):
        if grid[row][x] == num:
            return False
            
    # Check column
    for x in range(8):
        if grid[x][col] == num:
            return False
            
    # Check horizontal constraints
    h_cons = [(1,6), (3,6), (4,7), (7,4)]
    for r,c in h_cons:
        if row == r:
            if col == c and grid[r][c+1] != 0:
                if num >= grid[r][c+1]:
                    return False
            if col == c+1 and grid[r][c] != 0:
                if num <= grid[r][c]:
                    return False
                    
    # Check vertical constraints
    v_cons = [(0,1,'v'), (1,1,'^'), (1,7,'^'), (2,7,'v'), (4,4,'^'), (5,3,'v')]
    for r,c,t in v_cons:
        if col == c:
            if row == r and grid[r+1][c] != 0:
                if t == 'v' and num >= grid[r+1][c]:
                    return False
                if t == '^' and num <= grid[r+1][c]:
                    return False
            if row == r+1 and grid[r][c] != 0:
                if t == 'v' and grid[r][c] >= num:
                    return False
                if t == '^' and grid[r][c] <= num:
                    return False
                    
    return True

def find_empty(grid):
    # Find cell with minimum possible values first
    min_possibilities = 9
    best_pos = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                count = sum(1 for num in range(1,9) if is_valid(grid, i, j, num))
                if count < min_possibilities:
                    min_possibilities = count
                    best_pos = (i,j)
                    if count == 1:  # Can't get better than this
                        return best_pos
    
    return best_pos

def solve(grid):
    pos = find_empty(grid)
    if not pos:
        return True
        
    row, col = pos
    for num in range(1,9):
        if is_valid(grid, row, col, num):
            grid[row][col] = num
            if solve(grid):
                return True
            grid[row][col] = 0
            
    return False

# Initialize grid
grid = [[0]*8 for _ in range(8)]

# Set initial values
initials = [
    (0,1,2), (0,2,1), (0,4,5), (0,5,7),
    (1,1,1), (1,5,3),
    (2,0,4), (2,2,5),
    (3,1,7), (3,2,4), (3,3,6), (3,7,5),
    (4,1,3), (4,2,7), (4,6,1),
    (5,3,2), (5,6,5),
    (6,0,6), (6,4,2), (6,5,4),
    (7,4,4), (7,5,1), (7,6,7)
]

for r,c,v in initials:
    grid[r][c] = v

if solve(grid):
    print_puzzle(grid)
else:
    print("No solution exists")