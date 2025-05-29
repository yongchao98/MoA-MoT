from collections import deque

# Define the grid
grid = [
    ">>>>>>>>",
    ">eeee>e>",
    ">ee>>>>>",
    ">eeeeee>",
    ">e>ee>>e",
    ">>ez>3e>",
    ">eee>e>e",
    ">eeeee>e",
    ">>>>>>>>"
]

# Define the start and goal positions
start = (5, 6)
goal = (4, 2)

# Define the possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue for BFS
queue = deque([start])
# Initialize the visited set
visited = set([start])
# Initialize the steps dictionary
steps = {start: 0}

# Perform BFS
while queue:
    current = queue.popleft()
    if current == goal:
        break
    for move in moves:
        next_pos = (current[0] + move[0], current[1] + move[1])
        if (0 <= next_pos[0] < len(grid) and 0 <= next_pos[1] < len(grid[0]) and
                grid[next_pos[0]][next_pos[1]] == 'e' and next_pos not in visited):
            queue.append(next_pos)
            visited.add(next_pos)
            steps[next_pos] = steps[current] + 1

# Output the number of steps to reach the goal
print(steps[goal])