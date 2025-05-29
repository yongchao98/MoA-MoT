from collections import deque

# Define the grid
grid = [
    "bbbbbbbbb",
    "bpppppppb",
    "bpbjppbpb",
    "bpppbpbpb",
    "bpbpbpppb",
    "bbbpppppb",
    "bpbbCbbbb",
    "bpppbpppb",
    "bbbbbbbbb"
]

# Find the start position 'j'
start = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == 'j':
            start = (i, j)
            break
    if start:
        break

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS to find the shortest path
def bfs(start):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # Check if we reached the goal 'C'
        if grid[x][y] == 'C':
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] in ('p', 'C'):  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path found

# Find the minimum steps
min_steps = bfs(start)
print(min_steps)