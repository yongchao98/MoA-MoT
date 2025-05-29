from collections import deque

# Define the grid
grid = [
    "$$$$$$$$",
    "$$rrrrr$",
    "$rrrr$$$",
    "$$/rr$r$",
    "$rrr$r$$",
    "$$rrr$$$",
    "$$rrCr$$",
    "$$$$$$$$"
]

# Define the start and goal positions
start = (3, 2)  # Row 3, Column 2
goal = (6, 4)   # Row 6, Column 4

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# BFS initialization
queue = deque([(start, 0)])  # (position, steps)
visited = set()
visited.add(start)

# BFS loop
while queue:
    (current_row, current_col), steps = queue.popleft()
    
    # Check if we reached the goal
    if (current_row, current_col) == goal:
        print(steps)
        break
    
    # Explore neighbors
    for move in moves:
        new_row, new_col = current_row + move[0], current_col + move[1]
        if (0 <= new_row < len(grid) and 0 <= new_col < len(grid[0]) and
            grid[new_row][new_col] != '$' and (new_row, new_col) not in visited):
            visited.add((new_row, new_col))
            queue.append(((new_row, new_col), steps + 1))