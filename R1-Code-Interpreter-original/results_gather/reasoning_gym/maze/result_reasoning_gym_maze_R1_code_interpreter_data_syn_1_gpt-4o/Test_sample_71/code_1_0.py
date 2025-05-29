from collections import deque

# Define the grid
grid = [
    "||||||||||",
    "||W||W}WW|",
    "|WWW|W|WW|",
    "|W||WWWW||",
    "|WWW`WW|||",
    "|||W|WWWW|",
    "||WW|WW|W|",
    "||WWW||WW|",
    "|WWWW|W|W|",
    "||||||||||"
]

# Define the start and goal positions
start = (1, 6)  # Position of '}'
goal = (4, 4)   # Position of '`'

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue with the start position and step count
queue = deque([(start, 0)])
visited = set([start])

# Perform BFS
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we have reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore all possible moves
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] == 'W' and (nx, ny) not in visited:
            queue.append(((nx, ny), steps + 1))
            visited.add((nx, ny))