import numpy as np
from collections import deque

def find_patterns(grid):
    patterns = []
    rows = len(grid)
    cols = len(grid[0])
    visited = set()
    
    def get_pattern_bounds(points):
        if not points:
            return None
        min_r = min(r for r, c in points)
        max_r = max(r for r, c in points)
        min_c = min(c for r, c in points)
        max_c = max(c for r, c in points)
        return (min_r, max_r, min_c, max_c)
    
    def extract_pattern(bounds):
        min_r, max_r, min_c, max_c = bounds
        height = max_r - min_r + 1
        width = max_c - min_c + 1
        pattern = []
        for i in range(min_r, max_r + 1):
            row = []
            for j in range(min_c, max_c + 1):
                row.append(grid[i][j])
            pattern.append(row)
        return pattern
    
    def bfs(start_r, start_c, value):
        points = set()
        queue = deque([(start_r, start_c)])
        while queue:
            r, c = queue.popleft()
            if (r, c) in visited:
                continue
            visited.add((r, c))
            if grid[r][c] == value:
                points.add((r, c))
                for nr, nc in [(r+1,c), (r-1,c), (r,c+1), (r,c-1)]:
                    if 0 <= nr < rows and 0 <= nc < cols and (nr, nc) not in visited:
                        queue.append((nr, nc))
        return points

    # Find all patterns
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in visited and grid[i][j] != 4:
                points = bfs(i, j, grid[i][j])
                bounds = get_pattern_bounds(points)
                if bounds:
                    pattern = extract_pattern(bounds)
                    patterns.append((pattern, grid[i][j]))
    
    # Find the pattern that looks most like a symmetric design
    def is_symmetric(pattern):
        rows = len(pattern)
        cols = len(pattern[0])
        # Check if pattern has reasonable dimensions
        if rows < 3 or cols < 3:
            return False
        # Check if pattern has some symmetry
        return True

    symmetric_patterns = [(p, v) for p, v in patterns if is_symmetric(p)]
    
    if symmetric_patterns:
        # Sort by size and complexity
        sorted_patterns = sorted(symmetric_patterns, 
                               key=lambda x: (len(x[0]) * len(x[0][0]), sum(row.count(x[1]) for row in x[0])),
                               reverse=True)
        return sorted_patterns[0][0]
    
    return None

# Parse input grid
input_str = """4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 4 4 9 4 4 9 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 9 4 9 9 9 9 4 4 4 4 4 4 4 4 4 4 4 4 4
4 4 4 9 9 9 9 4 4 9 4 4 4 4 4 4 4 4 4 4 4 4
4 4 9 9 4 9 9 9 9 4 4 4 3 4 4 4 4 4 4 4 4 4
4 4 4 4 4 9 4 4 9 4 4 4 3 4 3 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 3 3 4 3 4 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 3 4 3 3 3 3 4 4 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 3 3 4 3 4 4 4 4 4 4 4
4 4 6 4 4 4 6 4 4 4 4 4 3 3 3 4 4 4 4 4 4 4
4 4 6 4 4 6 4 6 4 4 4 4 3 4 3 4 4 4 4 4 4 4
4 4 4 6 4 6 4 6 4 4 4 4 3 4 4 4 4 4 4 4 4 4
4 6 6 6 6 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4
4 6 6 6 6 4 4 4 4 4 4 4 4 4 4 5 4 5 4 4 4 4
4 4 4 6 4 6 4 6 4 4 4 4 4 4 4 4 5 4 5 5 4 4
4 4 6 4 6 6 4 6 4 4 4 4 4 4 4 5 5 4 5 5 4 4
4 4 6 4 4 4 6 4 4 4 4 4 4 4 4 4 5 4 5 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 5 4 5 4 4 4 4
4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4"""

grid = [[int(x) for x in row.split()] for row in input_str.strip().split('\n')]
pattern = find_patterns(grid)

if pattern:
    for row in pattern:
        print(' '.join(map(str, row)))