from collections import deque

# Define the grid
grid = [
    "yyyyyyyy",
    "yy666yyy",
    "yyy666}y",
    "y666666y",
    "y666666y",
    "yy6666yy",
    "y66y6$yy",
    "yyyyyyyy"
]

# Define the start and goal positions
start = (2, 6)  # Position of '}'
goal = (6, 5)   # Position of '$'

# Define possible movements (up, down, left, right)
movements = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue with the start position and step count
queue = deque([(start, 0)])
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
    for dx, dy in movements:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] == '6' or grid[nx][ny] == '$':
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))