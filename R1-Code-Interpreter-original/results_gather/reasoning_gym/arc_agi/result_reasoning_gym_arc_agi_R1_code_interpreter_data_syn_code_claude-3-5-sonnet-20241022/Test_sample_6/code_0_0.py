def find_patterns(input_grid, output_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    
    # Function to check if a position is part of a valid pattern
    def is_pattern(r, c):
        # Check if we can form a 3x3 pattern around this position
        if r < 1 or r >= rows-1 or c < 1 or c >= cols-1:
            return False
        
        # Check if this position is a 9
        if input_grid[r][c] != 9:
            return False
            
        # Check surrounding positions for 6s
        pattern = []
        for i in range(r-1, r+2):
            for j in range(c-1, c+2):
                if input_grid[i][j] not in [6, 9]:
                    return False
                pattern.append(input_grid[i][j])
        return True
    
    # Find all patterns and their preservation in output
    patterns = []
    for i in range(rows):
        for j in range(cols):
            if is_pattern(i, j):
                preserved = all(input_grid[i+di][j+dj] == output_grid[i+di][j+dj] 
                              for di in [-1,0,1] for dj in [-1,0,1])
                patterns.append((i, j, preserved))
    
    return patterns

# Test input grids
test_input = [
    [8,8,8,8,8,8,8,8,6,6,6,8],
    [8,6,6,6,6,6,6,8,6,6,6,8],
    [8,6,9,6,6,6,6,8,8,8,8,8],
    [8,6,6,6,6,9,6,8,6,6,6,6],
    [8,6,6,6,6,6,6,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,6,9,6],
    [8,8,8,8,6,6,6,8,6,6,6,6],
    [8,8,8,8,6,6,9,8,8,8,8,8],
    [8,8,8,8,6,6,6,8,8,8,8,8],
    [6,6,6,8,9,6,6,8,8,8,8,8],
    [9,6,6,8,8,8,8,8,6,9,6,9],
    [6,6,6,8,8,8,8,8,6,6,6,6],
    [8,8,8,8,8,8,8,8,6,9,6,6]
]

# Print the patterns found
print("Analyzing patterns...")
print("\nInput grid patterns:")
for row in test_input:
    print(' '.join(map(str, row)))