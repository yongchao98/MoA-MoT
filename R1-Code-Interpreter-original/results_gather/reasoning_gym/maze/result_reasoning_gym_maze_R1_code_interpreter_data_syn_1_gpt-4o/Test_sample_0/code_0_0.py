from collections import deque

# Define the grid
grid = [
    "999999999",
    "9Cr9rrrr9",
    "9rrrr9r99",
    "9rr9rirr9",
    "9rrrrr9r9",
    "9rrrrrr99",
    "9r99rr999",
    "9rrr9r9r9",
    "999999999"
]

# Find the start position 'i'
start = None
for r in range(len(grid)):
    for c in range(len(grid[0])):
        if grid[r][c] == 'i':
            start = (r, c)
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
        (r, c), steps = queue.popleft()
        
        # Check if we reached the goal
        if grid[r][c] == 'C':
            return steps
        
        # Explore neighbors
        for dr, dc in directions:
            nr, nc = r + dr, c + dc
            if 0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and (nr, nc) not in visited:
                if grid[nr][nc] in ('r', 'C'):  # Can move to passage or goal
                    visited.add((nr, nc))
                    queue.append(((nr, nc), steps + 1))
    
    return -1  # If goal is not reachable

# Find the minimum number of steps
min_steps = bfs(start)
print(min_steps)