from collections import deque

# Define the grid
grid = [
    "ppppppppp",
    "p``p`p``p",
    "p`p`p`p`p",
    "pp`pp```p",
    "pp`p````p",
    "p```p`;pp",
    "p````p``p",
    "p``d```pp",
    "ppppppppp"
]

# Convert grid to a list of lists for easier manipulation
grid = [list(row) for row in grid]

# Define start and goal positions
start = (5, 6)
goal = (7, 3)

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
                if grid[nx][ny] == '`' or (nx, ny) == goal:
                    visited.add((nx, ny))
                    queue.append(((nx, ny), steps + 1))
    
    return -1  # If no path is found

# Find the minimum number of steps
min_steps = bfs(start, goal)
print(min_steps)