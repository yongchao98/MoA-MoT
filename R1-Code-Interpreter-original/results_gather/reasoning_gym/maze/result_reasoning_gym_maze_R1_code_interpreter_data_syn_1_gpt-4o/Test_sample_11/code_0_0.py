from collections import deque

# Define the grid
grid = [
    "4444444",
    "4AAAAi4",
    "4A4A4A4",
    "4A4AA44",
    "44AAAA4",
    "44A`444",
    "4444444"
]

# Define the start and goal positions
start = (5, 3)  # Position of '`'
goal = (1, 5)   # Position of 'i'

# Define possible movements: up, down, left, right
movements = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue for BFS and visited set
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# Perform BFS
while queue:
    (current_x, current_y), steps = queue.popleft()
    
    # Check if we have reached the goal
    if (current_x, current_y) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in movements:
        new_x, new_y = current_x + dx, current_y + dy
        if (0 <= new_x < len(grid) and 0 <= new_y < len(grid[0]) and
            grid[new_x][new_y] in 'Ai' and (new_x, new_y) not in visited):
            visited.add((new_x, new_y))
            queue.append(((new_x, new_y), steps + 1))