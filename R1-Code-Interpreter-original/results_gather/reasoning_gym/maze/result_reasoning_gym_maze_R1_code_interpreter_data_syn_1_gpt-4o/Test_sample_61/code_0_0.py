from collections import deque

# Define the grid
grid = [
    '""""""""',
    '"++"++"',
    '++"+++"',
    '"++++J"',
    '"++++"+',
    '++}+++"',
    '""+"+++"',
    '""""""""'
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Define the start and goal positions
start = (5, 2)  # Position of '}'
goal = (3, 5)   # Position of 'J'

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize BFS
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# Perform BFS
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we have reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        
        # Check if the new position is within bounds and not a wall
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and grid[nx][ny] != '"' and (nx, ny) not in visited:
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))