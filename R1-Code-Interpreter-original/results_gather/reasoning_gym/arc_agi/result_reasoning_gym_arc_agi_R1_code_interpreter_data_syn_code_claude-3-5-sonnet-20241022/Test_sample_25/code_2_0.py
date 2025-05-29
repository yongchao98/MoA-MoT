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
    # Convert input string to 2D list of integers
    grid = []
    for line in input_grid.split('\n'):
        row = [int(x) for x in line.strip().split()]
        grid.append(row)
    grid = np.array(grid)
    output = np.copy(grid)
    
    # Find all clusters
    clusters = find_clusters(grid)
    
    # Process the grid according to the rules
    for i in range(10):
        for j in range(10):
            # Handle 0s - merge and move to corners
            if grid[i][j] == 0:
                output[i][j] = 1  # Clear original positions
            # Keep 5s and 6s in place
            elif grid[i][j] in [5, 6]:
                output[i][j] = grid[i][j]
            # Maintain structure of 2s
            elif grid[i][j] == 2:
                output[i][j] = 2
            else:
                output[i][j] = 1
    
    # Place merged clusters
    # Move 0s to bottom-right
    zero_count = np.sum(grid == 0)
    if zero_count > 0:
        for i in range(min(3, zero_count)):
            for j in range(min(3, zero_count)):
                output[9-i][9-j] = 0
    
    # Keep 6s in their relative position
    six_positions = np.where(grid == 6)
    if len(six_positions[0]) > 0:
        for i, j in zip(six_positions[0], six_positions[1]):
            output[i][j] = 6
            
    # Keep 5s in their relative position
    five_positions = np.where(grid == 5)
    if len(five_positions[0]) > 0:
        for i, j in zip(five_positions[0], five_positions[1]):
            output[i][j] = 5
    
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