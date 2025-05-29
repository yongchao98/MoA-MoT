def print_grid(grid):
    for row in grid:
        print(','.join(row))

def is_valid(grid, r, c, val, size=7):
    # Check row
    if val in grid[r]:
        return False
    
    # Check column
    if val in [grid[i][c] for i in range(size) if grid[i][c] != '.']:
        return False
    
    # Check minor diagonal if this position is on it
    if r + c == size - 1:
        for i in range(size):
            if grid[i][size-1-i] != '.' and grid[i][size-1-i] != val:
                return False
    
    return True

def solve(grid, r=0, c=0, size=7):
    if c == size:
        r += 1
        c = 0
    if r == size:
        return True
    
    # Skip if cell is already filled
    if grid[r][c] != '.':
        return solve(grid, r, c+1)
    
    # If on minor diagonal, try only the letter that should be there
    if r + c == size - 1:
        # Find the letter that should be on minor diagonal
        diag_letter = None
        for i in range(size):
            if grid[i][size-1-i] != '.':
                diag_letter = grid[i][size-1-i]
                break
        
        if diag_letter is not None:
            if is_valid(grid, r, c, diag_letter):
                grid[r][c] = diag_letter
                if solve(grid, r, c+1):
                    return True
            grid[r][c] = '.'
            return False
    
    # Try each possible letter
    for val in 'abcdefg':
        if is_valid(grid, r, c, val):
            grid[r][c] = val
            if solve(grid, r, c+1):
                return True
    grid[r][c] = '.'
    return False

# Initialize grid
initial_grid = [
    ['f', 'd', 'a', '.', '.', 'g', '.'],
    ['d', '.', '.', 'b', '.', '.', 'f'],
    ['.', 'c', 'b', '.', 'e', '.', '.'],
    ['.', '.', 'g', 'e', '.', '.', '.'],
    ['.', 'g', 'e', 'f', '.', 'a', 'c'],
    ['.', 'e', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'a', '.', '.', '.']
]

if solve(initial_grid):
    print_grid(initial_grid)
else:
    print("No solution exists")