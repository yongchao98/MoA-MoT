def find_connected_groups(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    groups = []
    
    def dfs(r, c, group):
        if (r, c) in visited or r < 0 or r >= rows or c < 0 or c >= cols or grid[r][c] != 8:
            return
        visited.add((r, c))
        group.append((r, c))
        for dr, dc in [(0,1), (1,0), (0,-1), (-1,0)]:
            dfs(r + dr, c + dc, group)
    
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 8 and (i,j) not in visited:
                group = []
                dfs(i, j, group)
                if group:
                    groups.append(group)
    return groups

def should_convert_to_five(group, grid):
    # Get group dimensions
    min_r = min(r for r,c in group)
    max_r = max(r for r,c in group)
    min_c = min(c for r,c in group)
    max_c = max(c for r,c in group)
    height = max_r - min_r + 1
    width = max_c - min_c + 1
    
    # If group is too large or too small, don't convert
    if height > 5 or width > 5 or len(group) < 3:
        return False
    
    # Check if group forms a pattern (roughly rectangular or symbol-like)
    pattern_density = len(group) / (height * width)
    return pattern_density >= 0.4 and height >= 2 and width >= 2

def process_grid(input_grid):
    grid = [row[:] for row in input_grid]
    groups = find_connected_groups(grid)
    
    for group in groups:
        if should_convert_to_five(group):
            for r, c in group:
                grid[r][c] = 5
    
    return grid

# Parse input grid
input_str = """4 4 4 4 4 4 8 8 8 8 8 4 4 4 4 4
8 8 8 8 4 4 8 4 4 8 4 4 4 4 4 4
4 4 4 8 4 4 4 4 4 8 4 4 4 4 4 4
4 4 4 4 4 4 4 4 8 8 8 8 8 8 8 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
8 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
8 8 8 8 4 4 4 4 4 4 4 4 4 4 4 4
8 4 4 8 4 4 4 4 8 8 8 8 8 8 4 4
8 4 4 8 4 4 4 4 4 4 8 4 4 8 4 4
8 4 4 8 4 4 4 4 4 4 8 4 4 8 4 4
8 4 4 8 4 4 4 4 4 4 8 8 8 8 4 4
8 8 8 8 4 4 8 8 8 4 4 4 4 8 4 4
4 4 4 4 4 4 4 4 8 4 4 4 4 4 4 4
4 4 4 4 4 8 4 4 8 4 4 4 4 4 4 8
4 4 4 4 4 8 8 8 8 4 4 4 4 4 4 8"""

input_grid = [[int(x) for x in line.split()] for line in input_str.split('\n')]
output_grid = process_grid(input_grid)

# Print result
for row in output_grid:
    print(' '.join(map(str, row)))