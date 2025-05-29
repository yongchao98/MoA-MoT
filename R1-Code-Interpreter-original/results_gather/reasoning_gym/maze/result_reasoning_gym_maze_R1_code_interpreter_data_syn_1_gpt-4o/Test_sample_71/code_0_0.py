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

# Define start and goal positions
start = (1, 7)
goal = (4, 3)

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue for BFS
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# BFS algorithm
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we have reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore all possible moves
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        
        # Check if the new position is within bounds and not a wall
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] != '|' and (nx, ny) not in visited:
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))