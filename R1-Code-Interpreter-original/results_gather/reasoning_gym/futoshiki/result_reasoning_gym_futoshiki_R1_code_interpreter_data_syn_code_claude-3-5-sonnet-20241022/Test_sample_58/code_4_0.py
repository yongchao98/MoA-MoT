def print_solution(grid):
    for i in range(5):
        row = []
        for j in range(5):
            row.append(str(grid[i][j]))
        print("   ".join(row))

def get_possible_values(grid, row, col):
    if grid[row][col] != 0:
        return set()
    
    possible = set(range(1, 6))
    
    # Remove values from same row
    for j in range(5):
        if grid[row][j] != 0:
            possible.discard(grid[row][j])
    
    # Remove values from same column
    for i in range(5):
        if grid[i][col] != 0:
            possible.discard(grid[i][col])
    
    # Apply vertical constraints
    if col == 1:  # Column 2
        if row == 1:  # Row 2
            if grid[2][1] != 0:  # Value below
                possible = {x for x in possible if x > grid[2][1]}
        if row == 2:  # Row 3
            if grid[1][1] != 0:  # Value above
                possible = {x for x in possible if x < grid[1][1]}
    
    if col == 4:  # Column 5
        if row == 0:  # Row 1
            if grid[1][4] != 0:  # Value below
                possible = {x for x in possible if x > grid[1][4]}
        if row == 1:  # Row 2
            if grid[0][4] != 0:  # Value above
                possible = {x for x in possible if x < grid[0][4]}
    
    return possible

def find_best_cell(grid):
    min_len = 6
    best_cell = None
    best_values = None
    
    for i in range(5):
        for j in range(5):
            if grid[i][j] == 0:
                possible = get_possible_values(grid, i, j)
                if len(possible) < min_len:
                    min_len = len(possible)
                    best_cell = (i, j)
                    best_values = possible
                if len(possible) == 0:
                    return None, None  # No valid values for this cell
    
    return best_cell, best_values

def solve(grid):
    cell, values = find_best_cell(grid)
    
    if cell is None and values is None:
        # Check if grid is complete
        for i in range(5):
            for j in range(5):
                if grid[i][j] == 0:
                    return False
        return True
    
    if not values:  # No valid values for selected cell
        return False
    
    row, col = cell
    for val in sorted(values):  # Try values in ascending order
        grid[row][col] = val
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

# Validate initial grid
valid = True
for i in range(5):
    for j in range(5):
        if grid[i][j] != 0:
            temp = grid[i][j]
            grid[i][j] = 0
            possible = get_possible_values(grid, i, j)
            if temp not in possible:
                valid = False
                print(f"Invalid initial value {temp} at position ({i+1}, {j+1})")
            grid[i][j] = temp

if valid and solve(grid):
    print_solution(grid)
else:
    print("No solution exists")