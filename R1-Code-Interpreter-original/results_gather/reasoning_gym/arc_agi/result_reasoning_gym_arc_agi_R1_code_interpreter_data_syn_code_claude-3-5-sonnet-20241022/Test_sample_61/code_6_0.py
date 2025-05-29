import numpy as np
from collections import deque

def get_connected_component(grid, start_x, start_y, value):
    rows, cols = len(grid), len(grid[0])
    component = set()
    queue = deque([(start_x, start_y)])
    
    while queue:
        x, y = queue.popleft()
        if (x, y) in component:
            continue
        if x < 0 or x >= rows or y < 0 or y >= cols:
            continue
        if grid[x][y] != value:
            continue
            
        component.add((x, y))
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            queue.append((x + dx, y + dy))
    
    return component

def expand_region(grid, output, start_points, value, blocked):
    rows, cols = len(grid), len(grid[0])
    queue = deque(list(start_points))
    visited = set(start_points)
    
    while queue:
        x, y = queue.popleft()
        output[x][y] = value
        
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = x + dx, y + dy
            if ((nx, ny) not in visited and (nx, ny) not in blocked and
                0 <= nx < rows and 0 <= ny < cols and
                (grid[nx][ny] == 7 or grid[nx][ny] == 2)):
                visited.add((nx, ny))
                queue.append((nx, ny))
    
    return visited

def solve_grid(grid):
    rows, cols = len(grid), len(grid[0])
    output = np.full((rows, cols), 9)  # Default background
    processed = set()
    
    # First process special numbers (non-2 and non-7)
    special_regions = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] not in [2, 7]:
                special_regions.append((grid[i][j], {(i, j)}))
    
    # Sort special regions by value (to ensure consistent processing order)
    special_regions.sort(key=lambda x: x[0])
    
    # Process 2's regions
    two_components = []
    for i in range(rows):
        for j in range(cols):
            if grid[i][j] == 2 and (i, j) not in processed:
                component = get_connected_component(grid, i, j, 2)
                processed.update(component)
                
                # Find nearest special number
                min_dist = float('inf')
                nearest_value = 0
                
                for x, y in component:
                    for value, region in special_regions:
                        for sx, sy in region:
                            dist = abs(x - sx) + abs(y - sy)
                            if dist < min_dist:
                                min_dist = dist
                                nearest_value = value
                
                two_components.append((nearest_value, component))
    
    # Sort two_components by value
    two_components.sort(key=lambda x: x[0])
    
    # Combine all regions
    all_regions = special_regions + two_components
    
    # Process each region in order
    blocked = set()
    for value, region in all_regions:
        expanded = expand_region(grid, output, region, value, blocked)
        blocked.update(expanded)
    
    return output

# Test input
test_input = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 0, 7, 6, 7],
    [7, 7, 2, 2, 2, 2, 2, 7, 7, 2, 7, 7, 7],
    [7, 7, 2, 7, 7, 7, 2, 7, 7, 2, 2, 7, 7],
    [7, 2, 2, 7, 7, 2, 2, 7, 7, 7, 2, 2, 2],
    [7, 2, 7, 7, 7, 2, 7, 7, 7, 7, 5, 7, 7],
    [7, 2, 7, 9, 7, 2, 7, 7, 7, 7, 7, 7, 7],
    [7, 2, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7],
    [7, 8, 7, 7, 7, 2, 2, 2, 7, 7, 7, 7, 7],
    [7, 2, 7, 7, 7, 7, 7, 2, 7, 7, 7, 7, 7],
    [7, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
    [7, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 2, 2, 7, 7, 7, 3, 7, 7, 7],
    [7, 7, 7, 7, 7, 2, 7, 7, 7, 7, 7, 7, 7]
]

result = solve_grid(test_input)
print('\n'.join(' '.join(map(str, row)) for row in result))