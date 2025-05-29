def find_connected_regions(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    regions = []
    
    def get_region(r, c, val, distance=1):
        points = set()
        for i in range(max(0, r-distance), min(rows, r+distance+1)):
            for j in range(max(0, c-distance), min(cols, c+distance+1)):
                if grid[i][j] == val and (i,j) not in visited:
                    points.add((i,j))
                    visited.add((i,j))
                    new_points = get_region(i, j, val, distance)
                    points.update(new_points)
        return points

    # First get 8s (corners)
    corner_points = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 8:
                corner_points.append((i,j))
                visited.add((i,j))
    
    # Then get other numbers
    for val in [7, 4, 2, 1]:  # Process in priority order
        for i in range(rows):
            for j in range(cols):
                if grid[i][j] == val and (i,j) not in visited:
                    region = get_region(i, j, val, distance=2)  # Increased merge distance
                    if region:
                        regions.append((val, region))
    
    return corner_points, regions

def create_output(input_grid):
    corners, regions = find_connected_regions(input_grid)
    
    # Determine output size based on content
    out_height = 11  # Based on example outputs
    out_width = 11
    
    output = [[0 for _ in range(out_width)] for _ in range(out_height)]
    
    # Place corners
    if corners:
        output[0][0] = 8
        output[0][-1] = 8
        output[-1][0] = 8
        output[-1][-1] = 8
    
    # Process each region
    for val, points in regions:
        min_r = min(r for r,c in points)
        max_r = max(r for r,c in points)
        min_c = min(c for r,c in points)
        max_c = max(c for r,c in points)
        
        # Calculate block dimensions
        height = max_r - min_r + 1
        width = max_c - min_c + 1
        
        # Calculate position in output grid
        if min_r < len(input_grid) // 2:  # Top half
            out_r = 2
        else:  # Bottom half
            out_r = out_height - height - 2
            
        if min_c < len(input_grid[0]) // 2:  # Left half
            out_c = 2
        else:  # Right half
            out_c = out_width - width - 2
        
        # Place the block
        for i in range(height):
            for j in range(width):
                if 0 <= out_r + i < out_height and 0 <= out_c + j < out_width:
                    output[out_r + i][out_c + j] = val
    
    # Print result
    result = ""
    for row in output:
        result += " ".join(map(str, row)) + "\n"
    print(result.strip())

# Test input
test_input = [
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,7,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,7,7,7,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,4,0,2,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,4,4,4,4,0,0,0,0,2,2,2,2,0,0],
    [0,0,0,0,4,4,4,4,0,0,0,0,2,2,2,2,0,0],
    [0,0,0,0,4,4,4,4,0,0,0,0,2,2,2,2,0,0],
    [0,0,0,0,4,4,4,4,0,0,0,0,2,2,2,2,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8]
]

create_output(test_input)