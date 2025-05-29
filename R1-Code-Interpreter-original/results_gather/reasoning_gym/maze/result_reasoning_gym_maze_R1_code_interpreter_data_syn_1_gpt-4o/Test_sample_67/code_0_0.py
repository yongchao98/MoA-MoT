from collections import deque

# Define the grid
grid = [
    ['P', 'P', 'P', 'P', 'P', 'P'],
    ['P', 'P', 'M', 'P', 'M', 'P'],
    ['P', 'M', 'M', 'M', 'M', 'P'],
    ['P', 'M', 'M', 'P', 'p', 'P'],
    ['P', 'M', '/', 'M', 'P', 'P'],
    ['P', 'P', 'P', 'P', 'P', 'P']
]

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Find the starting point 'p'
start = None
for i in range(len(grid)):
    for j in range(len(grid[0])):
        if grid[i][j] == 'p':
            start = (i, j)
            break
    if start:
        break

# BFS to find the shortest path
def bfs(start):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)

    while queue:
        (x, y), steps = queue.popleft()

        # Check if we have reached the goal
        if grid[x][y] == '/':
            return steps

        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] in ('M', '/'):
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))

    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start)
print(min_steps)