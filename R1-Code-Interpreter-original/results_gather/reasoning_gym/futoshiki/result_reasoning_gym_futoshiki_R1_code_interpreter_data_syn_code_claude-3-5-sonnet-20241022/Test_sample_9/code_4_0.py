def get_domain(grid, row, col, h_constraints, v_constraints):
    domain = set(range(1, 8))
    
    # Remove numbers already in row
    for j in range(7):
        if grid[row][j] != 0:
            domain.discard(grid[row][j])
    
    # Remove numbers already in column
    for i in range(7):
        if grid[i][col] != 0:
            domain.discard(grid[i][col])
    
    # Check horizontal constraints
    if col > 0 and h_constraints[row][col-1] == '<' and grid[row][col-1] != 0:
        domain = {x for x in domain if x > grid[row][col-1]}
    if col > 0 and h_constraints[row][col-1] == '>' and grid[row][col-1] != 0:
        domain = {x for x in domain if x < grid[row][col-1]}
    if col < 6 and h_constraints[row][col] == '<' and grid[row][col+1] != 0:
        domain = {x for x in domain if x < grid[row][col+1]}
    if col < 6 and h_constraints[row][col] == '>' and grid[row][col+1] != 0:
        domain = {x for x in domain if x > grid[row][col+1]}
    
    # Check vertical constraints
    if row > 0 and v_constraints[row-1][col] == '∧' and grid[row-1][col] != 0:
        domain = {x for x in domain if x > grid[row-1][col]}
    if row > 0 and v_constraints[row-1][col] == '∨' and grid[row-1][col] != 0:
        domain = {x for x in domain if x < grid[row-1][col]}
    if row < 6 and v_constraints[row][col] == '∧' and grid[row+1][col] != 0:
        domain = {x for x in domain if x < grid[row+1][col]}
    if row < 6 and v_constraints[row][col] == '∨' and grid[row+1][col] != 0:
        domain = {x for x in domain if x > grid[row+1][col]}
    
    return sorted(list(domain))

def get_next_cell(grid):
    min_domain_size = 8
    next_cell = None
    next_domain = None
    
    for i in range(7):
        for j in range(7):
            if grid[i][j] == 0:
                domain = get_domain(grid, i, j, h_constraints, v_constraints)
                if len(domain) < min_domain_size:
                    min_domain_size = len(domain)
                    next_cell = (i, j)
                    next_domain = domain
    
    return next_cell, next_domain

# Initialize puzzle
grid = [
    [0, 0, 5, 0, 4, 0, 0],
    [0, 0, 0, 0, 0, 6, 4],
    [7, 0, 0, 0, 0, 4, 0],
    [0, 2, 0, 0, 7, 1, 0],
    [0, 0, 0, 0, 3, 0, 0],
    [0, 1, 3, 2, 0, 0, 0],
    [6, 4, 0, 0, 1, 0, 0]
]

h_constraints = [
    ['', '', '', '<', '', ''],
    ['<', '', '', '<', '', ''],
    ['', '', '', '', '<', ''],
    ['', '', '', '', '', ''],
    ['', '', '>', '', '', ''],
    ['', '', '', '', '', ''],
    ['', '', '<', '>', '', '']
]

v_constraints = [
    ['∧', '', '∧', '∨', '∧', ''],
    ['∧', '', '', '∧', '', ''],
    ['', '', '∧', '', '', '∨'],
    ['', '', '', '', '', ''],
    ['', '', '', '', '', ''],
    ['∧', '', '', '', '', '']
]

# Solve using stack to avoid recursion
stack = []
while True:
    next_cell, domain = get_next_cell(grid)
    
    if not next_cell:  # Solution found
        break
        
    if not domain:  # No valid values for this cell
        if not stack:  # No previous states to backtrack to
            print("No solution exists")
            exit()
        grid, next_cell, domain = stack.pop()
        continue
        
    row, col = next_cell
    value = domain[0]
    domain = domain[1:]
    
    # Save current state for backtracking
    if domain:
        stack.append((
            [row[:] for row in grid],
            next_cell,
            domain
        ))
    
    grid[row][col] = value

# Print solution
for i in range(7):
    for j in range(7):
        print(grid[i][j], end='')
        if j < 6:
            print(' ' + (h_constraints[i][j] if h_constraints[i][j] else ' ') + ' ', end='')
    print()
    if i < 6:
        for j in range(7):
            print((v_constraints[i][j] if v_constraints[i][j] else ' ') + '   ', end='')
        print()