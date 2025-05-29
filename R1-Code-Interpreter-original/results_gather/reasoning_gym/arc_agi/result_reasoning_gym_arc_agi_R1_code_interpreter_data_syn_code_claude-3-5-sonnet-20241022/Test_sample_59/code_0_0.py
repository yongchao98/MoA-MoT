def find_pattern(grid):
    # Find non-zero regions
    n = len(grid)
    patterns = []
    visited = set()
    
    def get_region(i, j, visited):
        if (i, j) in visited or i < 0 or i >= n or j < 0 or j >= n or grid[i][j] == 0:
            return []
        region = [(i, j)]
        visited.add((i, j))
        for ni, nj in [(i+1,j), (i-1,j), (i,j+1), (i,j-1)]:
            region.extend(get_region(ni, nj, visited))
        return region
    
    # Find all connected regions
    for i in range(n):
        for j in range(n):
            if grid[i][j] != 0 and (i,j) not in visited:
                region = get_region(i, j, visited)
                if region:
                    patterns.append(region)
    return patterns

def transform_grid(input_grid):
    n = len(input_grid)
    output = [[0]*n for _ in range(n)]
    patterns = find_pattern(input_grid)
    
    for pattern in patterns:
        # Find bounds of pattern
        min_i = min(i for i,j in pattern)
        max_i = max(i for i,j in pattern)
        min_j = min(j for i,j in pattern)
        max_j = max(j for i,j in pattern)
        
        # Find frame and inner numbers
        frame_num = None
        inner_num = None
        for i, j in pattern:
            if i in (min_i, max_i) or j in (min_j, max_j):
                frame_num = input_grid[i][j]
            else:
                inner_num = input_grid[i][j]
        
        # Expand pattern
        new_min_i = max(0, min_i - 1)
        new_max_i = min(n-1, max_i + 1)
        new_min_j = max(0, min_j - 1)
        new_max_j = min(n-1, max_j + 1)
        
        # Fill expanded pattern
        for i in range(new_min_i, new_max_i + 1):
            for j in range(new_min_j, new_max_j + 1):
                if i in (new_min_i, new_max_i) or j in (new_min_j, new_max_j):
                    output[i][j] = frame_num
                else:
                    output[i][j] = inner_num if inner_num else frame_num

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