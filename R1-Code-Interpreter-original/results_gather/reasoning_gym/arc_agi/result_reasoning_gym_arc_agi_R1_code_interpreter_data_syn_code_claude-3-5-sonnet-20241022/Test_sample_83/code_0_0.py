def find_connected_regions(grid):
    if not grid or not grid[0]:
        return []
    
    rows, cols = len(grid), len(grid[0])
    visited = set()
    regions = []
    
    def dfs(r, c, value):
        if (r, c) in visited or r < 0 or r >= rows or c < 0 or c >= cols or grid[r][c] != value:
            return 0
        
        visited.add((r, c))
        size = 1
        
        # Check all 4 directions
        for dr, dc in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            size += dfs(r + dr, c + dc, value)
            
        return size
    
    # Find all regions
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in visited and grid[i][j] != 3:  # Skip background (3)
                size = dfs(i, j, grid[i][j])
                regions.append((size, grid[i][j]))
    
    return regions

# Parse input grid
input_str = """3 3 3 3 3 3 3 3 3 3 7 7 7 7 3
9 9 9 9 9 9 9 9 3 3 7 7 7 7 3
9 9 9 9 9 9 9 9 3 3 7 7 7 7 3
9 9 9 9 9 9 9 9 3 3 7 7 7 7 3
9 9 9 9 9 9 9 9 3 3 7 7 7 7 3
9 9 9 9 9 9 9 9 3 3 3 3 3 3 3
9 9 9 9 9 9 9 9 3 2 2 2 2 3 3
3 3 3 3 3 3 3 3 3 2 2 2 2 3 3
3 6 6 6 3 3 3 3 3 2 2 2 2 3 3
3 6 6 6 5 5 5 5 3 2 2 2 2 3 3
3 6 6 6 5 5 5 5 3 3 3 3 3 3 3
3 3 3 3 5 5 5 5 3 3 4 4 4 4 3
3 3 3 3 5 5 5 5 3 3 4 3 3 4 3
3 3 3 3 3 3 3 3 3 3 4 4 4 4 3"""

# Convert input string to grid
grid = [[int(x) for x in row.split()] for row in input_str.split('\n')]

# Find all regions and their sizes
regions = find_connected_regions(grid)

# Sort regions by size (descending) and value
regions.sort(reverse=True)

# Print the value of the largest region
print(regions[0][1])