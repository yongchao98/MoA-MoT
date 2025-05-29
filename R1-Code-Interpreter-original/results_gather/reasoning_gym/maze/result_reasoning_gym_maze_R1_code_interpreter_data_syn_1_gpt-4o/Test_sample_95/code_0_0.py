from collections import deque

# Define the grid
grid = [
    "5555555555",
    "5<5<<<5<55",
    "5<5<<<55<5",
    "5<<-55<<<5",
    "555<<<<<<5",
    "5<5<<<<<<5",
    "5<5<55<<55",
    "5<+<<<<5<5",
    "5555<<5<<5",
    "5555555555"
]

# Find the start and goal positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == '-':
            start = (i, j)
        elif cell == '+':
            goal = (i, j)

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS setup
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
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] in '<+':  # Can move to passage or goal
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))