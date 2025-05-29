from collections import deque

# Define the grid
grid = [
    "////////",
    "//((//(/",
    "/((((((/",
    "/(((/(}/",
    "/((/(((/",
    "/((((///",
    "/((((C//",
    "////////"
]

# Define the start and goal
start = (3, 7)  # Position of '}'
goal = (6, 5)   # Position of 'C'

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the BFS queue
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# Perform BFS
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] == '(' or grid[nx][ny] == 'C':
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))