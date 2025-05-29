from collections import deque

def find_shortest_path(grid):
    # Define the grid
    grid = [
        ['O', 'X', 'O', 'X', 'X', 'O', 'X'],
        ['O', 'O', 'X', 'X', 'X', 'O', '#'],
        ['X', 'O', 'X', 'O', 'X', 'O', 'O'],
        ['O', 'O', 'X', 'X', 'O', 'X', 'X'],
        ['O', 'O', 'X', 'X', 'O', 'O', 'O'],
        ['O', 'O', 'O', 'O', '*', 'X', 'O']
    ]

    # Directions for movement: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    direction_names = ['up', 'down', 'left', 'right']

    # Find the start and destination points
    start = None
    destination = None
    for i in range(len(grid)):
        for j in range(len(grid[0])):
            if grid[i][j] == '*':
                start = (i, j)
            elif grid[i][j] == '#':
                destination = (i, j)

    # Initialize the BFS queue
    queue = deque([(start, [])])
    visited = set()
    visited.add(start)

    # Perform BFS
    while queue:
        (current, path) = queue.popleft()

        # Check if we reached the destination
        if current == destination:
            return ' '.join(path)

        # Explore neighbors
        for d, name in zip(directions, direction_names):
            ni, nj = current[0] + d[0], current[1] + d[1]
            if 0 <= ni < len(grid) and 0 <= nj < len(grid[0]) and grid[ni][nj] in ('O', '#') and (ni, nj) not in visited:
                visited.add((ni, nj))
                queue.append(((ni, nj), path + [name]))

    return "infeasible"

# Find and print the shortest path
print(find_shortest_path([]))