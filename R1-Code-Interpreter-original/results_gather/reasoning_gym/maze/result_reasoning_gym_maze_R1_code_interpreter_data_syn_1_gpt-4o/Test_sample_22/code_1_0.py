from collections import deque

# Define the grid
grid = [
    "hhhhhhhh",
    "hhUUUhhh",
    "hhUUUh*h",
    "hUUUUUUh",
    "hUhhUUUh",
    "hCUhUhhh",
    "hUhUUUUh",
    "hhhhhhhh"
]

# Define start and goal positions
start = (2, 6)
goal = (5, 1)

# Define possible moves: up, down, left, right
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue for BFS
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# Perform BFS
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we have reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and
            grid[nx][ny] == 'U' and (nx, ny) not in visited):
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))