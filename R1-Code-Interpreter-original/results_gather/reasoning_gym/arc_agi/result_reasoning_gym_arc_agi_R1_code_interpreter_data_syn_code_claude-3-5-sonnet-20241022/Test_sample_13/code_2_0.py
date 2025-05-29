def find_connected_groups(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    groups = []
    
    def dfs(r, c, group):
        if (r, c) in visited or r < 0 or r >= rows or c < 0 or c >= cols or grid[r][c] != 8:
            return
        visited.add((r, c))
        group.append((r, c))
        for dr, dc in [(0,1), (1,0), (0,-1), (-1,0), (1,1), (1,-1), (-1,1), (-1,-1)]:
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
    min_r = min(r for r,c in group)
    max_r = max(r for r,c in group)
    min_c = min(c for r,c in group)
    max_c = max(c for r,c in group)
    height = max_r - min_r + 1
    width = max_c - min_c + 1
    
    # Size constraints
    if len(group) < 3 or len(group) > 25:
        return False
    if height > 7 or width > 7:
        return False
    
    # Check if the pattern is isolated
    for r, c in group:
        for dr, dc in [(0,1), (1,0), (0,-1), (-1,0)]:
            nr, nc = r + dr, c + dc
            if 0 <= nr < len(grid) and 0 <= nc < len(grid[0]):
                if grid[nr][nc] == 8 and (nr, nc) not in group:
                    return False
    
    # Check pattern characteristics
    pattern_density = len(group) / (height * width)
    return (pattern_density >= 0.3 and 
            height >= 2 and width >= 2 and 
            (height >= 3 or width >= 3))

def process_grid(input_grid):
    grid = [row[:] for row in input_grid]
    groups = find_connected_groups(grid)
    
    # Sort groups by size (process larger groups first)
    groups.sort(key=len, reverse=True)
    
    for group in groups:
        if should_convert_to_five(group, grid):
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

for row in output_grid:
    print(' '.join(map(str, row)))