from collections import deque

# Define the grid
grid = [
    "ssssssssss",
    "saaaWsaaas",
    "ssassssaas",
    "ssaaaassss",
    "saaaaassss",
    "saaaaa@aas",
    "sasaaasass",
    "saaaaaaaas",
    "saaasaasss",
    "ssssssssss"
]

# Define start and goal positions
start = (5, 6)
goal = (1, 4)

# Define possible movements: up, down, left, right
movements = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS initialization
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# BFS loop
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in movements:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and (nx, ny) not in visited and grid[nx][ny] in 'aW':
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))