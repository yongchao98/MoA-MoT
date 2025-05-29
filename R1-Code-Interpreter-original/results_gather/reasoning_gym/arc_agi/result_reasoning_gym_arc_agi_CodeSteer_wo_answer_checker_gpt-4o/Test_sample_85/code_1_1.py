from collections import deque

def transform_grid(input_grid):
    # Find the positions of '1' and '9'
    pos_1 = None
    pos_9 = None
    for i, row in enumerate(input_grid):
        for j, value in enumerate(row):
            if value == 1:
                pos_1 = (i, j)
            elif value == 9:
                pos_9 = (i, j)
    
    # Create a copy of the input grid for the output
    output_grid = [row[:] for row in input_grid]
    
    # Directions for moving in the grid (right, left, down, up)
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]
    
    # BFS to find the path from '1' to '9'
    if pos_1 and pos_9:
        queue = deque([pos_1])
        visited = set([pos_1])
        parent = {pos_1: None}
        
        while queue:
            current = queue.popleft()
            if current == pos_9:
                break
            for direction in directions:
                neighbor = (current[0] + direction[0], current[1] + direction[1])
                if (0 <= neighbor[0] < len(input_grid) and
                    0 <= neighbor[1] < len(input_grid[0]) and
                    neighbor not in visited):
                    queue.append(neighbor)
                    visited.add(neighbor)
                    parent[neighbor] = current
        
        # Reconstruct the path from '1' to '9' and fill it with '8'
        step = pos_9
        while step != pos_1:
            step = parent[step]
            if step != pos_1:
                output_grid[step[0]][step[1]] = 8
    
    return output_grid

# Test input grid
input_grid = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 9, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
]

output_grid = transform_grid(input_grid)

# Print the output grid
for row in output_grid:
    print(' '.join(map(str, row)))