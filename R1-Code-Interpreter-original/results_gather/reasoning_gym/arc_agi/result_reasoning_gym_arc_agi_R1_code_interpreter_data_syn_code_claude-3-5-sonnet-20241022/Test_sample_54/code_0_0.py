def find_patterns(grid):
    patterns = []
    rows = len(grid)
    cols = len(grid[0])
    
    # Look for 2x2 or larger patterns of non-7 numbers
    for i in range(rows-1):
        for j in range(cols-1):
            if grid[i][j] != 7:
                # Find extent of pattern
                height = 1
                width = 1
                while i + height < rows and grid[i + height][j] != 7:
                    height += 1
                while j + width < cols and grid[i][j + width] != 7:
                    width += 1
                if height > 1 or width > 1:
                    patterns.append({
                        'pos': (i, j),
                        'value': grid[i][j],
                        'height': height,
                        'width': width
                    })
    return patterns

def generate_output(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output = [[7 for _ in range(cols)] for _ in range(rows)]
    
    # Find patterns in input
    patterns = find_patterns(input_grid)
    
    # For each pattern found, apply transformation rules
    for pattern in patterns:
        i, j = pattern['pos']
        val = pattern['value']
        h = pattern['height']
        w = pattern['width']
        
        # If pattern is 2x2 or larger of 9s
        if val == 9 and (h >= 2 or w >= 2):
            # Place vertical bars of 3s
            for r in range(4):
                output[i+r][j-3] = 3
                output[i+r][j+w+2] = 3
            # Place horizontal bars of 9s
            for c in range(3):
                output[i-2][j-1+c] = 9
                output[i-2][j+w+c] = 9
                output[i+h+1][j-1+c] = 9
                output[i+h+1][j+w+c] = 9
        
        # If pattern is vertical line of 3s
        elif val == 3 and h > 1:
            # Create symmetric pattern
            center = j
            for offset in [-4, 4]:
                for r in range(3):
                    output[i+r][center+offset] = 3
            for c in range(3):
                output[i-2][center-1+c] = 9
                output[i+h+1][center-1+c] = 9
    
    return output

# Test input
test_input = [
    [7]*16 for _ in range(10)] + \
    [[7,7,7,7,7,7,7,9,9,7,7,7,7,7,7,7]] + \
    [[7]*16 for _ in range(5)]

# Generate and print output
output = generate_output(test_input)
print("Output grid:")
for row in output:
    print(" ".join(map(str, row)))