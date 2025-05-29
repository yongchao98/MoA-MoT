def find_7_regions(grid):
    rows = len(grid)
    cols = len(grid[0])
    visited = set()
    regions = []
    
    def dfs(r, c):
        if (r, c) not in visited and 0 <= r < rows and 0 <= c < cols and grid[r][c] == 7:
            visited.add((r, c))
            region = [(r, c)]
            for dr, dc in [(0,1), (1,0), (0,-1), (-1,0), (1,1), (-1,-1), (1,-1), (-1,1)]:
                nr, nc = r + dr, c + dc
                region.extend(dfs(nr, nc) or [])
            return region
        return []

    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 7 and (i, j) not in visited:
                region = dfs(i, j)
                if region:
                    regions.append(region)
    return regions

def generate_output(input_grid):
    rows = len(input_grid)
    cols = len(input_grid[0])
    output = [[2 for _ in range(cols)] for _ in range(rows)]
    
    # Copy 7s from input
    for i in range(rows):
        for j in range(cols):
            if input_grid[i][j] == 7:
                output[i][j] = 7
    
    # Find 7 regions
    regions = find_7_regions(input_grid)
    
    # Add shadow effect (0s)
    for region in regions:
        # Find region boundaries
        min_r = min(r for r, _ in region)
        max_r = max(r for r, _ in region)
        min_c = min(c for _, c in region)
        max_c = max(c for _, c in region)
        
        # Add top shadow
        if min_r > 0:
            for c in range(min_c, max_c + 1):
                if min_r - 1 >= 0:
                    output[min_r - 1][c] = 0
        
        # Add diagonal shadows
        for r, c in region:
            # Right diagonal shadow
            if c < cols - 1 and r > 0 and input_grid[r][c+1] != 7:
                output[r-1][c+1] = 0
            # Left diagonal shadow
            if c > 0 and r > 0 and input_grid[r][c-1] != 7:
                output[r-1][c-1] = 0

    return output

def print_grid(grid):
    result = []
    for row in grid:
        result.append(' '.join(str(x) for x in row))
    print('\n'.join(result))

# Test input
input_grid = [
    [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
    [2,2,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2,2,2],
    [2,7,7,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2,2],
    [7,7,7,7,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2],
    [7,7,7,7,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2],
    [7,7,7,7,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2],
    [7,7,7,7,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2],
    [7,7,7,7,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2],
    [2,7,7,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2,2],
    [2,2,7,7,7,7,7,7,7,2,2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,2,7,7,7,7,2,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,7,7,7,7,7,7,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,7,7,7,7,7,7,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,7,7,7,7,7,7,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,2,7,7,7,7,2,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
    [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
]

output = generate_output(input_grid)
print_grid(output)