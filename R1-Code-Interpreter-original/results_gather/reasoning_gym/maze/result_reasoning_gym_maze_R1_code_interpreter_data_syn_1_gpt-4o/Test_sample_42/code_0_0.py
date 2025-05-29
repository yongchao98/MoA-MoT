from collections import deque

# Define the grid
grid = [
    ";;;;;;;;;;",
    ";;;$$;$$$;",
    ";;$$$$;$$;",
    ";$$$$$$$$;",
    ";$;E;$$;;;",
    ";$$;;$$$$;",
    ";$$;$$;$;;",
    ";$$;$;;$$;",
    ";;$$$;$e$;",
    ";;;;;;;;;;"
]

# Define start and goal positions
start = (8, 7)
goal = (4, 2)

# Define possible moves (up, down, left, right)
moves = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Initialize the queue for BFS
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
    
    # Explore all possible moves
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(grid) and 0 <= ny < len(grid[0]) and (nx, ny) not in visited:
            if grid[nx][ny] == '$' or grid[nx][ny] == 'E':
                visited.add((nx, ny))
                queue.append(((nx, ny), steps + 1))