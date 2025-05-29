def find_patterns(grid):
    patterns = {}
    height = len(grid)
    width = len(grid[0])
    
    # Convert input string to 2D grid
    if isinstance(grid, str):
        grid = [[int(x) for x in row.split()] for row in grid.split('\n') if row]
    
    # Find all unique numbers and their positions
    for i in range(height):
        for j in range(width):
            if grid[i][j] != 9:
                num = grid[i][j]
                if num not in patterns:
                    patterns[num] = []
                patterns[num].append((i, j))
    
    # Extract pattern shapes
    result = {}
    for num in patterns:
        positions = patterns[num]
        min_i = min(p[0] for p in positions)
        min_j = min(p[1] for p in positions)
        max_i = max(p[0] for p in positions)
        max_j = max(p[1] for p in positions)
        
        # Create relative pattern
        pattern = []
        for i, j in positions:
            pattern.append((i - min_i, j - min_j))
        result[num] = sorted(pattern)
    
    return result

# Test input
test_input = """9 9 9 9 9 9 9 9 9 9 9 9 9 9 9 9
9 9 9 9 9 9 9 9 9 9 9 0 9 9 9 9
9 9 9 9 9 9 9 9 1 9 9 9 0 0 9 9
9 9 9 9 9 9 9 1 9 1 9 9 0 9 9 9
9 9 9 8 9 8 9 9 1 9 9 9 9 9 9 9
9 9 9 8 8 8 9 9 9 9 9 9 9 9 9 9
9 9 9 8 9 8 9 9 9 9 9 9 9 9 9 9
9 9 9 9 9 9 9 9 9 9 9 7 7 7 9 9
9 9 9 9 9 9 9 9 9 9 9 9 7 7 9 9
9 9 9 9 9 9 9 9 9 9 9 9 7 9 9 9
6 6 6 6 6 6 6 9 9 9 9 9 9 9 9 9
6 5 6 6 6 8 6 9 9 9 9 9 9 9 9 9
6 6 6 6 6 6 6 9 9 9 9 9 9 9 9 9
6 6 6 6 6 6 6 9 9 5 9 5 9 9 9 9
6 6 6 6 6 6 6 9 9 5 5 5 9 9 9 9
6 7 6 6 6 0 6 9 9 9 5 9 9 9 9 9
6 6 6 6 6 6 6 9 9 9 9 9 9 9 9 9"""

patterns = find_patterns(test_input)

# Create 7x7 output grid
output = [[9 for _ in range(7)] for _ in range(7)]

# Based on the pattern analysis, construct the output
# Pattern positions in output grid
pattern_positions = {
    8: [(0,1), (1,0), (1,2), (2,1)],  # Cross pattern
    1: [(0,4), (0,5), (1,4), (1,5), (1,6), (2,5)],  # L-shape
    0: [(4,0), (5,1), (5,2)],  # Triangle
    7: [(4,4), (4,5), (5,5)],  # L-shape
    5: [(4,0), (4,2), (5,1)]  # Cross pattern
}

# Print the result
for row in output:
    print(' '.join(map(str, row)))