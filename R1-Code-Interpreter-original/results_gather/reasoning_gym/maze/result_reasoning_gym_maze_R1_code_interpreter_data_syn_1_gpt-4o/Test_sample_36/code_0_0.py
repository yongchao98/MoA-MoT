from collections import deque

# Define the grid
grid = [
    "XXXXXXXX",
    "XXX___eX",
    "X_X_XX_X",
    "X____XXX",
    "X_X__#_X",
    "XX_X___X",
    "X___X_XX",
    "XXXXXXXX"
]

# Define the directions for movement: right, left, down, up
directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

# Find the start and goal positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == '#':
            start = (i, j)
        elif cell == 'e':
            goal = (i, j)

# BFS to find the shortest path
def bfs(start, goal):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)

    while queue:
        (x, y), steps = queue.popleft()

        # Check if we have reached the goal
        if (x, y) == goal:
            return steps

        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] != 'X' and (nx, ny) not in visited:
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))

    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)