from collections import deque

# Define the grid
grid = [
    "hhhhhhhh",
    "hhUUUhhh",
    "hhUUUh*h",
    "hUUUUUUh",
    "hUhhUUUh",
    "hCUhUhhh",
    "hUhUUUUh",
    "hhhhhhhh"
]

# Define the start and goal positions
start = (2, 6)  # Row 2, Column 6
goal = (5, 1)   # Row 5, Column 1

# Define the possible movements (up, down, left, right)
movements = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue for BFS and the visited set
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# Perform BFS
while queue:
    (current_row, current_col), steps = queue.popleft()
    
    # Check if we have reached the goal
    if (current_row, current_col) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for dr, dc in movements:
        new_row, new_col = current_row + dr, current_col + dc
        if (0 <= new_row < len(grid) and 0 <= new_col < len(grid[0]) and
            grid[new_row][new_col] != 'h' and (new_row, new_col) not in visited):
            visited.add((new_row, new_col))
            queue.append(((new_row, new_col), steps + 1))