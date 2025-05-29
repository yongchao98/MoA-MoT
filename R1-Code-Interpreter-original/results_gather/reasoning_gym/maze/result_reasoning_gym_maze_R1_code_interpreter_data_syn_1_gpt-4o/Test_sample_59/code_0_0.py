from collections import deque

# Define the grid
grid = [
    "aaaaaaaaaa",
    "a.aaa.a..a",
    "a..a.aa.aa",
    "a.m..a.aaa",
    "aa.....aaa",
    "a...a....a",
    "aa..a.aa.a",
    "a^a..aaa.a",
    "a....aaa.a",
    "aaaaaaaaaa"
]

# Define the start and goal positions
start = (3, 2)  # 'm' is at (3, 2)
goal = (7, 1)   # '^' is at (7, 1)

# Define the possible movements (up, down, left, right)
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

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
    
    # Explore neighbors
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and (nx, ny) not in visited and grid[nx][ny] != 'a':
            visited.add((nx, ny))
            queue.append(((nx, ny), steps + 1))