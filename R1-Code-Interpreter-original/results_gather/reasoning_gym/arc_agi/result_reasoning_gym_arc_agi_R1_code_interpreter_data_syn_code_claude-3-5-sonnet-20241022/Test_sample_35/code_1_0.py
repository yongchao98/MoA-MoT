def should_transform_to_five(grid, i, j):
    rows, cols = len(grid), len(grid[0])
    
    # Check if current position is 0
    if grid[i][j] != '0':
        return False
    
    # Helper function to safely get grid value
    def get_val(r, c):
        if 0 <= r < rows and 0 <= c < cols:
            return grid[r][c]
        return None
    
    # Check 2x2 square patterns
    square_patterns = [
        [(0,0), (0,1), (1,0), (1,1)],
        [(-1,-1), (-1,0), (0,-1), (0,0)],
        [(-1,0), (-1,1), (0,0), (0,1)],
        [(0,-1), (0,0), (1,-1), (1,0)],
        [(0,0), (0,1), (1,0), (1,1)]
    ]
    
    for pattern in square_patterns:
        zeros = 0
        for dr, dc in pattern:
            val = get_val(i+dr, j+dc)
            if val == '0':
                zeros += 1
        if zeros >= 3:
            return True
    
    # Check cross pattern
    cross = [(i-1,j), (i+1,j), (i,j-1), (i,j+1)]
    cross_zeros = sum(1 for r,c in cross if get_val(r,c) == '0')
    if cross_zeros >= 2:
        return True
    
    return False

def transform_grid(input_grid):
    # Convert input string to 2D grid
    grid = [row.split() for row in input_grid.strip().split('\n')]
    output = [row[:] for row in grid]
    
    # Process each cell
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if should_transform_to_five(grid, i, j):
                output[i][j] = '5'
    
    return '\n'.join(' '.join(row) for row in output)

# Test input
test_input = """1 7 1 1 1 1 1 1 7 7 1 1 7 7 7 1 1 1
1 1 1 1 1 7 1 1 7 1 1 1 7 0 0 7 0 7
1 7 1 1 7 1 1 7 7 7 1 1 0 7 7 0 0 1
1 7 1 7 7 7 7 7 7 1 1 1 1 1 1 1 1 1
7 7 7 7 1 1 7 1 7 7 1 7 1 1 1 7 7 1
1 1 1 1 7 1 1 7 0 0 7 0 7 1 7 7 7 7
1 1 1 1 7 7 1 7 0 0 0 0 7 1 1 1 1 1
7 1 7 7 7 1 7 7 1 7 1 1 7 1 7 7 1 7
1 1 1 1 7 1 7 7 7 7 7 1 1 1 7 1 1 1
7 1 1 1 1 7 7 1 1 7 7 1 7 1 7 1 7 7
7 7 7 7 1 1 7 1 1 7 0 7 1 7 1 7 7 1
7 1 7 7 7 1 1 1 7 7 7 0 7 1 1 7 7 1
1 7 7 7 1 1 1 7 1 1 7 7 7 1 7 7 1 7
1 7 0 0 7 0 0 7 1 1 1 1 7 7 7 7 7 7
7 0 7 7 0 7 0 1 1 7 1 7 7 7 7 7 7 1
1 7 0 0 7 0 0 1 7 1 7 7 1 1 7 1 7 1
1 7 7 7 1 7 1 1 1 1 1 1 1 1 1 1 1 7"""

print(transform_grid(test_input))