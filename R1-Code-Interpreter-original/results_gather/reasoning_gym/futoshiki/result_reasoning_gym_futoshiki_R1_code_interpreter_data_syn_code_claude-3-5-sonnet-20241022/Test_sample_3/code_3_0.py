def print_solution(grid):
    h_constraints = [(0,2,'<'), (1,2,'>'), (1,5,'>'), (2,0,'>'), (3,1,'>'), (3,5,'<'), (5,2,'>'), (5,7,'>')]
    v_constraints = [(1,0,'∧'), (1,3,'∧'), (2,2,'∧'), (2,4,'∧'), (2,5,'∧'), (4,7,'∧'), (6,1,'∧'), (6,6,'∧')]
    
    result = ""
    for i in range(8):
        for j in range(8):
            result += str(grid[i][j]) + "   "
            if j < 7:
                found = False
                for r, c, symbol in h_constraints:
                    if r == i and c == j:
                        result += symbol + " "
                        found = True
                        break
                if not found:
                    result += "  "
        result += "\n"
        if i < 7:
            for j in range(8):
                found = False
                for r, c, symbol in v_constraints:
                    if r == i and c == j:
                        result += symbol + "   "
                        found = True
                        break
                if not found:
                    result += "    "
            result += "\n"
    print("<<<")
    print(result.rstrip())
    print(">>>")

def is_valid(grid, row, col, num):
    # Check row and column uniqueness
    for x in range(8):
        if grid[row][x] == num or grid[x][col] == num:
            return False
    
    # Horizontal constraints (<, >)
    h_pairs = [
        (0, 2, '<'), (1, 2, '>'), (1, 5, '>'),
        (2, 0, '>'), (3, 1, '>'), (3, 5, '<'),
        (5, 2, '>'), (5, 7, '>')
    ]
    
    # Vertical constraints (∧)
    v_pairs = [
        (1, 0), (1, 3), (2, 2), (2, 4), (2, 5),
        (4, 7), (6, 1), (6, 6)
    ]
    
    # Check horizontal constraints
    for r, c, symbol in h_pairs:
        if r == row:
            if c == col and symbol == '<':
                if col < 7 and grid[r][col+1] != 0 and num >= grid[r][col+1]:
                    return False
            elif c == col and symbol == '>':
                if col < 7 and grid[r][col+1] != 0 and num <= grid[r][col+1]:
                    return False
            elif c == col-1 and symbol == '<':
                if grid[r][c] != 0 and grid[r][c] >= num:
                    return False
            elif c == col-1 and symbol == '>':
                if grid[r][c] != 0 and grid[r][c] <= num:
                    return False

    # Check vertical constraints (∧ means upper number must be smaller)
    for r, c in v_pairs:
        if c == col:
            if row == r:  # Current cell is upper cell
                if row < 7 and grid[row+1][col] != 0 and num >= grid[row+1][col]:
                    return False
            elif row == r+1:  # Current cell is lower cell
                if grid[r][col] != 0 and grid[r][col] >= num:
                    return False
    
    return True

def find_best_empty(grid):
    min_possibilities = 9
    best_cell = None
    
    for i in range(8):
        for j in range(8):
            if grid[i][j] == 0:
                count = sum(1 for num in range(1, 9) if is_valid(grid, i, j, num))
                if count < min_possibilities:
                    min_possibilities = count
                    best_cell = (i, j)
    
    return best_cell

def solve(grid):
    empty = find_best_empty(grid)
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

# Initial grid
grid = [
    [1, 4, 5, 0, 0, 0, 0, 3],
    [2, 0, 0, 0, 7, 0, 0, 0],
    [4, 0, 0, 0, 0, 6, 7, 0],
    [0, 0, 0, 4, 0, 1, 0, 8],
    [0, 0, 6, 2, 0, 3, 0, 0],
    [0, 0, 0, 0, 3, 0, 0, 0],
    [0, 0, 0, 0, 8, 0, 0, 5],
    [3, 6, 0, 0, 4, 7, 0, 0]
]

if solve(grid):
    print_solution(grid)
else:
    print("No solution exists")