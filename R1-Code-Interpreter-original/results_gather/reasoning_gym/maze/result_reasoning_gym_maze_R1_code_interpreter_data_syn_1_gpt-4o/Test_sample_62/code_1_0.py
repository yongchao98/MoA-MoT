from collections import deque

# Define the grid
grid = [
    "(((((((((",
    "(n(((nnn(",
    "(nn(nnn((",
    "(nnnn(nn(",
    "(f(nnn(((",
    "(nnn(nn((",
    "((nnnnnn(",
    "((nnnnnz(",
    "(((((((((("
]

# Define start and goal positions
start = (4, 0)
goal = (7, 7)

# Define possible moves: up, down, left, right
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the BFS queue
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# Perform BFS
while queue:
    (x, y), steps = queue.popleft()
    
    # Check if we reached the goal
    if (x, y) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and (nx, ny) not in visited and grid[nx][ny] == 'n':
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))