def transform_grid(input_grid):
    n = len(input_grid)
    output = [[0]*n for _ in range(n)]
    
    def find_patterns():
        patterns = []
        visited = set()
        
        def get_pattern(i, j):
            if input_grid[i][j] == 0:
                return None
            
            # Find boundaries
            top = bottom = i
            left = right = j
            
            # Expand while non-zero
            while top > 0 and input_grid[top-1][j] != 0:
                top -= 1
            while bottom < n-1 and input_grid[bottom+1][j] != 0:
                bottom += 1
            while left > 0 and input_grid[i][left-1] != 0:
                left -= 1
            while right < n-1 and input_grid[i][right+1] != 0:
                right += 1
                
            return (top, left, bottom, right)
        
        for i in range(n):
            for j in range(n):
                if input_grid[i][j] != 0 and (i,j) not in visited:
                    pattern = get_pattern(i, j)
                    if pattern:
                        patterns.append(pattern)
                        for r in range(pattern[0], pattern[2]+1):
                            for c in range(pattern[1], pattern[3]+1):
                                visited.add((r,c))
        return patterns
    
    def expand_pattern(top, left, bottom, right):
        # Get frame and inner numbers
        frame_num = input_grid[top][left]
        inner_num = None
        for i in range(top+1, bottom):
            for j in range(left+1, right):
                if input_grid[i][j] != frame_num:
                    inner_num = input_grid[i][j]
                    break
            if inner_num:
                break
        
        # New boundaries
        new_top = max(0, top-1)
        new_left = max(0, left-1)
        new_bottom = min(n-1, bottom+1)
        new_right = min(n-1, right+1)
        
        # Fill expanded pattern
        for i in range(new_top, new_bottom+1):
            for j in range(new_left, new_right+1):
                if i == new_top or i == new_bottom or j == new_left or j == new_right:
                    output[i][j] = frame_num
                elif i >= top and i <= bottom and j >= left and j <= right:
                    # Preserve original pattern
                    output[i][j] = input_grid[i][j]
                else:
                    # Fill new inner area
                    output[i][j] = inner_num if inner_num else frame_num
    
    patterns = find_patterns()
    for pattern in patterns:
        expand_pattern(*pattern)
    
    return output

# Test input
test_input = [
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,6,6,6,0,0,0,0,0,0,0,0],
    [0,6,4,6,0,0,0,0,0,0,0,0],
    [0,6,6,6,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,6,6,6,6,0,0],
    [0,0,0,0,0,0,6,4,4,6,0,0],
    [0,0,0,0,0,0,6,4,4,6,0,0],
    [0,0,0,0,0,0,6,6,6,6,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0]
]

result = transform_grid(test_input)
for row in result:
    print(' '.join(map(str, row)))