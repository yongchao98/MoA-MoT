import numpy as np
from collections import defaultdict

def find_clusters(grid):
    rows, cols = len(grid), len(grid[0])
    visited = set()
    clusters = defaultdict(list)
    
    def dfs(r, c, num):
        if (r, c) in visited or r < 0 or r >= rows or c < 0 or c >= cols or grid[r][c] != num:
            return []
        visited.add((r, c))
        positions = [(r, c)]
        for dr, dc in [(0,1), (1,0), (0,-1), (-1,0)]:
            positions.extend(dfs(r+dr, c+dc, num))
        return positions
    
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in visited:
                cluster = dfs(i, j, grid[i][j])
                if cluster:
                    clusters[grid[i][j]].append(cluster)
    return clusters

def solve(input_grid):
    # Convert input grid to numpy array
    grid = np.array([list(map(int, row)) for row in input_grid.split('\n')])
    output = np.copy(grid)
    
    # Rules based on pattern analysis:
    # 1. Move clusters to corners or edges when possible
    # 2. Maintain cluster sizes
    # 3. Preserve relative positions of different numbers
    
    clusters = find_clusters(grid)
    
    # Process specific patterns
    for num in clusters:
        for cluster in clusters[num]:
            cluster_size = len(cluster)
            min_r = min(r for r, c in cluster)
            max_r = max(r for r, c in cluster)
            min_c = min(c for r, c in cluster)
            max_c = max(c for r, c in cluster)
            
            # If cluster is in top-left, move to top-right or bottom-left
            if min_r < 5 and min_c < 5:
                if num in [0, 8, 9]:  # Numbers that tend to move to corners
                    for r, c in cluster:
                        output[r][c] = grid[r][c]
                        if min_r <= 2 and min_c <= 2:
                            output[r][c] = grid[max_r-(r-min_r)][max_c-(c-min_c)]
            
            # If cluster is in middle, try to move to edges
            if 3 <= min_r <= 6 and 3 <= min_c <= 6:
                if num in [5, 6]:  # Numbers that tend to stay in place
                    continue
                for r, c in cluster:
                    if grid[r][c] == num:
                        output[r][c] = grid[r][c]
    
    # Format output
    result = '\n'.join(' '.join(str(x) for x in row) for row in output)
    print(result)

# Test input
test_input = """0 1 1 2 2 2 2 2 2 2
1 0 1 2 1 1 1 2 2 2
1 1 1 2 2 2 2 2 2 2
1 5 1 1 1 1 1 1 1 1
5 5 5 1 1 0 1 0 1 1
1 1 2 2 2 2 2 1 1 0
1 1 2 2 2 2 2 1 1 1
1 1 2 2 1 2 2 6 6 6
0 1 2 1 1 1 2 1 1 1
1 1 2 2 2 2 2 1 1 0"""

solve(test_input)