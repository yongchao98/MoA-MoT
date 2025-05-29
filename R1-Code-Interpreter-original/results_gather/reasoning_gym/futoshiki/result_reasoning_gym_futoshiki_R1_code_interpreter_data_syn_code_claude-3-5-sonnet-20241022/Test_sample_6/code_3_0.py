def check_constraints(grid, row, col, num, h_constraints, v_constraints):
    # Check row and column
    for i in range(9):
        if grid[row][i] == num or grid[i][col] == num:
            return False
            
    # Check horizontal constraints
    for r, c1, c2, type in h_constraints:
        if r == row:
            if c1 == col and grid[r][c2] != 0:
                if (type == '<' and num >= grid[r][c2]) or (type == '>' and num <= grid[r][c2]):
                    return False
            if c2 == col and grid[r][c1] != 0:
                if (type == '<' and grid[r][c1] >= num) or (type == '>' and grid[r][c1] <= num):
                    return False

    # Check vertical constraints
    for r1, r2, c, type in v_constraints:
        if c == col:
            if r1 == row and grid[r2][c] != 0:
                if (type == '∨' and num <= grid[r2][c]) or (type == '∧' and num >= grid[r2][c]):
                    return False
            if r2 == row and grid[r1][c] != 0:
                if (type == '∨' and grid[r1][c] <= num) or (type == '∧' and grid[r1][c] >= num):
                    return False
    return True

def find_empty(grid):
    min_possibilities = 10
    best_pos = None
    
    for i in range(9):
        for j in range(9):
            if grid[i][j] == 0:
                count = sum(1 for num in range(1, 10) 
                          if check_constraints(grid, i, j, num, h_constraints, v_constraints))
                if count < min_possibilities:
                    min_possibilities = count
                    best_pos = (i, j)
    return best_pos

def solve(grid, h_constraints, v_constraints):
    pos = find_empty(grid)
    if not pos:
        return True
    
    row, col = pos
    # Try numbers in an optimized order based on constraints
    numbers = list(range(1, 10))
    # Prioritize numbers based on constraints
    if any(r == row and c1 == col and type == '<' for r, c1, c2, type in h_constraints):
        numbers.sort()
    elif any(r == row and c1 == col and type == '>' for r, c1, c2, type in h_constraints):
        numbers.sort(reverse=True)
    
    for num in numbers:
        if check_constraints(grid, row, col, num, h_constraints, v_constraints):
            grid[row][col] = num
            if solve(grid, h_constraints, v_constraints):
                return True
            grid[row][col] = 0
    return False

# Initial grid
grid = [
    [8, 5, 9, 0, 4, 0, 1, 0, 0],
    [0, 4, 2, 5, 0, 0, 3, 0, 0],
    [7, 1, 0, 0, 0, 9, 0, 4, 0],
    [0, 0, 0, 0, 0, 6, 7, 3, 0],
    [3, 7, 0, 0, 9, 1, 0, 0, 0],
    [0, 0, 0, 2, 8, 0, 0, 0, 0],
    [6, 0, 7, 0, 0, 0, 0, 5, 0],
    [0, 9, 6, 0, 2, 0, 0, 0, 8],
    [0, 0, 5, 7, 0, 2, 0, 0, 0]
]

# Horizontal constraints
h_constraints = [
    (0, 2, 3, '<'), (0, 4, 5, '>'), (0, 5, 6, '>'), (0, 6, 7, '>'),
    (1, 0, 1, '<'), (1, 6, 7, '>'), (1, 7, 8, '>'),
    (2, 3, 4, '<'), (2, 5, 6, '>'), (2, 6, 7, '>'),
    (3, 1, 2, '>'), (3, 3, 4, '>'),
    (4, 2, 3, '>'), (4, 3, 4, '<'),
    (6, 4, 5, '<'), (6, 6, 7, '>'),
    (7, 2, 3, '>'), (7, 3, 4, '>'), (7, 6, 7, '>')
]

# Vertical constraints
v_constraints = [
    (1, 2, 4, '∨'),
    (3, 4, 3, '∨'),
    (3, 4, 6, '∨'),
    (4, 5, 6, '∧'),
    (6, 7, 0, '∨'),
    (7, 8, 1, '∨'),
    (3, 4, 8, '∧')
]

def print_solution(grid):
    for row in grid:
        print(' '.join(str(x) for x in row))

print("Solving Futoshiki puzzle...")
if solve(grid, h_constraints, v_constraints):
    print("Solution found:")
    print_solution(grid)
else:
    print("No solution exists")