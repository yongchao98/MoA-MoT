from collections import deque

# Define the grid
grid = [
    "999999999",
    "99//////9",
    "9/9////99",
    "9/9/99//9",
    "9/9///9/9",
    "9//////99",
    "9n9///I/9",
    "9/99///99",
    "999999999"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Find the start position 'n'
start = None
goal = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == 'n':
            start = (i, j)
        elif grid[i][j] == 'I':
            goal = (i, j)

# Directions for moving in the grid (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS to find the shortest path
def bfs(start, goal):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)

    while queue:
        (x, y), steps = queue.popleft()

        # Check if we reached the goal
        if (x, y) == goal:
            return steps

        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] in ('/', 'I'):  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))

# Get the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)