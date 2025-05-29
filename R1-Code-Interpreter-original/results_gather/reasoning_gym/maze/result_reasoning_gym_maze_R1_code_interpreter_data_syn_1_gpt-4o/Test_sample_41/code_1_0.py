from collections import deque

# Define the grid
grid = [
    "zzzzzz",
    "zffffz",
    "z#fzZz",
    "zfffzz",
    "zffzzz",
    "zzzzzz"
]

# Define start and goal positions
start = (2, 1)
goal = (2, 4)

# Define possible movements: up, down, left, right
movements = [(-1, 0), (1, 0), (0, -1), (0, 1)]

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
    for dx, dy in movements:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 6 and 0 <= ny < 6 and (nx, ny) not in visited and grid[nx][ny] == 'f':
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))