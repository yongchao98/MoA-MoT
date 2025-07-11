import collections

def find_safest_path():
    """
    Finds the least dangerous path from '@' to 'g' on the given map.

    The method uses Breadth-First Search (BFS) to find the shortest path,
    which in this specific map layout, is also the least dangerous as it moves
    away from the only major threat, the dragon 'D'. The path is then
    simplified as per the user's request.
    """
    # The map provided by the user.
    grid_str = """
/       - - - - - - -                                         
/       | . . . . |         ############                      
/       | . . . . |         #          #                      
/       | . g . . + ########### #        @                      
/       - - - - - - -         #        ---+- - -              
/                             #        | . . . . . |            
/                             #        | . ! . . . |            
/                             #        | . . . . . |            
/                             #        | . . . . . |            
/         - - - -             #        | . . . . . |            
/         | . . |             #######+ . . D . . |            
/         | < . + ###         #        | . . . . . |            
/         - - - -   #         #        | . ? . . . |            
/                   ####### #        - - - - - - -            
"""
    grid = [list(line) for line in grid_str.strip().split('\n')]
    rows = len(grid)
    cols = len(grid[0])

    # Find start ('@') and goal ('g') positions.
    start_pos, goal_pos = None, None
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '@':
                start_pos = (r, c)
            elif grid[r][c] == 'g':
                goal_pos = (r, c)

    if not start_pos or not goal_pos:
        print("Error: Start or goal not found on the map.")
        return

    # Define obstacles. The dragon 'D' is also treated as an impassable obstacle.
    obstacles = {'#', '-', '|', '/', 'D'}

    # Initialize a queue for BFS with (position, path_taken).
    queue = collections.deque([(start_pos, [])])
    visited = {start_pos}

    path_to_goal = None

    # Start BFS to find the shortest path.
    while queue:
        (r, c), path = queue.popleft()

        if (r, c) == goal_pos:
            path_to_goal = path
            break

        # Define possible moves: Up, Down, Left, Right.
        moves = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
        for move_char, (dr, dc) in moves.items():
            nr, nc = r + dr, c + dc

            # Check if the new position is valid.
            if 0 <= nr < rows and 0 <= nc < cols and \
               (nr, nc) not in visited and \
               grid[nr][nc] not in obstacles:
                
                visited.add((nr, nc))
                new_path = path + [move_char]
                queue.append(((nr, nc), new_path))
    
    if not path_to_goal:
        print("No path was found to the gold.")
        return

    # Simplify the path by removing consecutive duplicate moves.
    if not path_to_goal:
        simplified_path_str = ""
    else:
        simplified_path = [path_to_goal[0]]
        for i in range(1, len(path_to_goal)):
            if path_to_goal[i] != simplified_path[-1]:
                simplified_path.append(path_to_goal[i])
        simplified_path_str = "".join(simplified_path)
    
    print(simplified_path_str)

find_safest_path()