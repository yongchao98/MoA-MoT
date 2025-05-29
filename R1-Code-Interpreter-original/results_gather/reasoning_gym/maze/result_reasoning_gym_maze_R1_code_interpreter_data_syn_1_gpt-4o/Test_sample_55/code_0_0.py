from collections import deque

# Define the grid
grid = [
    "iiiiiiiii",
    "ioooooi(i",
    "iooiioioi",
    "ioooooooi",
    "iooioiioi",
    "ioooooioi",
    "ioioioidi",
    "ioiooiooi",
    "iiiiiiiii"
]

# Find the start position
start = None
for r, row in enumerate(grid):
    for c, char in enumerate(row):
        if char == '(':
            start = (r, c)
            break
    if start:
        break

# Define the goal character
goal = 'd'

# Define possible movements: up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS initialization
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# BFS loop
while queue:
    (r, c), steps = queue.popleft()
    
    # Check if we reached the goal
    if grid[r][c] == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dr, dc in directions:
        nr, nc = r + dr, c + dc
        if (0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and
            grid[nr][nc] != 'i' and (nr, nc) not in visited):
            visited.add((nr, nc))
            queue.append(((nr, nc), steps + 1))