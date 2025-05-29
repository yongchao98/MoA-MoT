from collections import deque

# Define the grid
grid = [
    "ZZZZZZZ",
    "Z-ZZ\\Z",
    "Z\\ZZ\\Z",
    "Z\\\\Z\\Z",
    "Z\\\\\\<Z",
    "Z\\\\Z\\Z",
    "ZZZZZZZ"
]

# Define the start and goal positions
start = (1, 1)  # Row 1, Column 1
goal = (4, 5)   # Row 4, Column 5

# Define possible moves: up, down, left, right
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue for BFS
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
    
    # Explore all possible moves
    for move in moves:
        new_row = current_row + move[0]
        new_col = current_col + move[1]
        new_position = (new_row, new_col)
        
        # Check if the new position is within bounds and is a passage
        if (0 <= new_row < len(grid) and 0 <= new_col < len(grid[0]) and
            grid[new_row][new_col] in ['\\', '<'] and new_position not in visited):
            visited.add(new_position)
            queue.append((new_position, steps + 1))