from collections import deque

# Define the grid
grid = [
    "vvvvvvvv",
    "vsssvssv",
    "vvsvssvv",
    "vsssvvsv",
    "vssssssv",
    "vssvsssv",
    "vssCv/vv",
    "vvvvvvvv"
]

# Define the start and goal positions
start = (6, 7)  # Position of '/'
goal = (6, 3)   # Position of 'C'

# Define possible movements: up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

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
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] == 's' or grid[nx][ny] == 'C':  # Can move to passage or goal
                queue.append(((nx, ny), steps + 1))
                visited.add((nx, ny))