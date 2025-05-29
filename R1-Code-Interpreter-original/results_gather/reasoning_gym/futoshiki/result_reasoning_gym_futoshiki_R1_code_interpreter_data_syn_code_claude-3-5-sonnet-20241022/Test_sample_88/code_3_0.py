def print_solution(grid):
    result = []
    # Define constraints positions
    h_less = [(1,6), (3,6), (4,7), (7,4)]  # positions of '<'
    v_down = [(0,1), (2,7), (5,3)]  # positions of '∨'
    v_up = [(1,1), (1,7), (4,4)]    # positions of '∧'
    
    for i in range(8):
        row = ""
        for j in range(7):
            row += str(grid[i][j])
            if (i,j) in h_less:
                row += " < "
            else:
                row += "   "
        row += str(grid[i][7])
        result.append(row)
        
        if i < 7:
            v_row = "    " * 8
            v_list = list(v_row)
            for r,c in v_down:
                if r == i:
                    v_list[c*4 + 1] = '∨'
            for r,c in v_up:
                if r == i:
                    v_list[c*4 + 1] = '∧'
            result.append("".join(v_list).rstrip())
    
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
    h_less = [(1,6), (3,6), (4,7), (7,4)]
    for r,c in h_less:
        if row == r:
            if col == c and grid[r][c+1] != 0:
                if num >= grid[r][c+1]:
                    return False
            if col == c+1 and grid[r][c] != 0:
                if num <= grid[r][c]:
                    return False
    
    # Check vertical constraints
    v_down = [(0,1), (2,7), (5,3)]
    v_up = [(1,1), (1,7), (4,4)]
    
    for r,c in v_down:
        if col == c:
            if row == r and grid[r+1][c] != 0:
                if num >= grid[r+1][c]:
                    return False
            if row == r+1 and grid[r][c] != 0:
                if grid[r][c] >= num:
                    return False
    
    for r,c in v_up:
        if col == c:
            if row == r and grid[r+1][c] != 0:
                if num <= grid[r+1][c]:
                    return False
            if row == r+1 and grid[r][c] != 0:
                if grid[r][c] <= num:
                    return False
    
    return True

def solve(grid):
    empty = None
    # Find first empty cell
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                empty = (i,j)
                break
        if empty:
            break
    
    if not empty:
        return True
    
    row, col = empty
    for num in range(1, 9):
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
    print_solution(grid)
else:
    print("No solution exists")