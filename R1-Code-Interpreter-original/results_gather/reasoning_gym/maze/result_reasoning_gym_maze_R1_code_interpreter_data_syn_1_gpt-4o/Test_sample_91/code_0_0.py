from collections import deque

# Define the grid
grid = [
    "6666666666",
    "6yyyyyyyy6",
    "66yyy66yy6",
    "66yyyyy6y6",
    "666yyyy6y6",
    "6yy6b6yQy6",
    "6y6yyyyyy6",
    "6y666yyyy6",
    "6yyy66yy66",
    "6666666666"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Directions for moving in the grid (right, left, down, up)
directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

# Find the starting point 'b'
start = None
for i in range(len(grid)):
    for j in range(len(grid[i])):
        if grid[i][j] == 'b':
            start = (i, j)
            break
    if start:
        break

# BFS to find the shortest path to 'Q'
def bfs(start):
    queue = deque([(start, 0)])  # (position, steps)
    visited = set()
    visited.add(start)
    
    while queue:
        (x, y), steps = queue.popleft()
        
        # Check if we've reached the goal
        if grid[x][y] == 'Q':
            return steps
        
        # Explore neighbors
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
                if grid[nx][ny] in ('y', 'Q'):  # Can move to passage or goal
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start)
print(min_steps)