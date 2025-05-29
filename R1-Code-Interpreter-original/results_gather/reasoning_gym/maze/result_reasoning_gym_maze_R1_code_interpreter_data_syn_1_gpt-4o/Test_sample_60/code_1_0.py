from collections import deque

# Define the grid
grid = [
    "RRRRRRRRRR",
    "RR%%RRR%%R",
    "RRRR%%%%%R",
    "R%%%%%RR%R",
    "R%%%%%%%%R",
    "R%%R%%%RRR",
    "R%%+%%%%%R",
    "R%%%R%%)RR",
    "R%%%%%RR%R",
    "RRRRRRRRRR"
]

# Define the start and goal positions
start = (6, 3)
goal = (7, 7)

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
    
    # Explore the neighbors
    for dx, dy in movements:
        nx, ny = x + dx, y + dy
        if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and
            grid[nx][ny] == '%' and (nx, ny) not in visited):
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))