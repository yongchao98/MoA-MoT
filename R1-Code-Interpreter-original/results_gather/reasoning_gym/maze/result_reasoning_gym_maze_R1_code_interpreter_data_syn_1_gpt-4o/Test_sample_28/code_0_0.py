from collections import deque

# Define the grid
grid = [
    "DDDDDDDDDD",
    "DDDAADAAAD",
    "DDAjAADDDD",
    "DADDADADAD",
    "DAAAAAAADD",
    "DAAAADAAAD",
    "DDADAADAAD",
    "DAAAADAAAD",
    "DAD,DDAAAD",
    "DDDDDDDDDD"
]

# Define start and goal positions
start = (8, 3)
goal = (2, 2)

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
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] == 'A' or (nx, ny) == goal:
                queue.append(((nx, ny), steps + 1))
                visited.add((nx, ny))