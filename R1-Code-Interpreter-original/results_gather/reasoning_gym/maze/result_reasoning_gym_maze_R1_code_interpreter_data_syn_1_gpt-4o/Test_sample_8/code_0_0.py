from collections import deque

# Define the grid
grid = [
    "WWWWWWWWW",
    "W--W--WWW",
    "W-W----WW",
    "W--4---WW",
    "W-------W",
    "W-W--WWWW",
    "W-WW--W-W",
    "WW-----&W",
    "WWWWWWWWW"
]

# Define start and goal positions
start = (7, 7)
goal = (3, 3)

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue with the start position and step count
queue = deque([(start, 0)])
visited = set()
visited.add(start)

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
        
        # Check if the new position is valid
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] != 'W' and (nx, ny) not in visited:
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))