from collections import deque

# Define the grid
grid = [
    "))))))",
    ")gGGG)",
    ")GGG))",
    ")G)GH)",
    ")G))G)",
    "))))))"
]

# Directions for moving in the grid (right, down, left, up)
directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]

# Find the start and goal positions
start = None
goal = None
for i, row in enumerate(grid):
    for j, cell in enumerate(row):
        if cell == 'g':
            start = (i, j)
        elif cell == 'H':
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
            if (0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and
                grid[nx][ny] != ')' and (nx, ny) not in visited):
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))

    return -1  # If no path is found

# Calculate the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)