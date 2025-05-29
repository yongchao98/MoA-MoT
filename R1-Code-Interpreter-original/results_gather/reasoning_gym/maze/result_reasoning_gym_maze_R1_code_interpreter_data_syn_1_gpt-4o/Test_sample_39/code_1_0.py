from collections import deque

# Define the grid
grid = [
    "FFFFFFFFFF",
    "F$$F$$F$$F",
    "FF$$$$$F$F",
    "FF$$F$$$FF",
    "F$$$FF$$FF",
    "F$$F$F$$<F",
    "F$F$$FF$$F",
    "F$$$F1$$$F",
    "FF$$F$FF$F",
    "FFFFFFFFFF"
]

# Define the start and goal positions
start = (5, 8)
goal = (7, 4)

# Define the possible movements (up, down, left, right)
movements = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue for BFS and the visited set
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
    for dx, dy in movements:
        nx, ny = x + dx, y + dy
        if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and
            grid[nx][ny] == '$' and (nx, ny) not in visited):
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))