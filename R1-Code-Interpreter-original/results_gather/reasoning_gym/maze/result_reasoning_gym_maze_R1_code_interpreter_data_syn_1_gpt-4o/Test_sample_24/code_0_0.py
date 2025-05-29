from collections import deque

# Define the grid
grid = [
    "*********",
    "*11**111*",
    "*@111111*",
    "*1*11111*",
    "*11***11*",
    "*11*1111*",
    "**11Y111*",
    "*1111111*",
    "*********"
]

# Define the start and goal positions
start = (6, 4)  # 'Y' position
goal = (2, 1)   # '@' position

# Define possible movements: up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue with the start position and step count
queue = deque([(start, 0)])
visited = set()
visited.add(start)

# BFS algorithm
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we have reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and
            grid[nx][ny] != '*' and (nx, ny) not in visited):
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))