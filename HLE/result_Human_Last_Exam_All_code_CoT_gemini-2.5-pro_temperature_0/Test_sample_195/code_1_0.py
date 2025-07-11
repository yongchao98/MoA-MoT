import collections

def solve_maze():
    """
    Solves the maze by finding the least dangerous path from '@' to 'g'.
    The least dangerous path is defined as the shortest path that avoids the dragon's room.
    """
    grid_str = [
        "/                                                                                ",
        "/  | . . . . . |          # # # # # # # # # # # #                                  ",
        "/  | . . . . . |          #                     #                                  ",
        "/  | . g . . + # # # # # # # #   @                                                ",
        "/  | . . . . . |            #   - - - + - - -                                      ",
        "/  - - - - - - -            #   | . . . . . . |                                    ",
        "/                           #   | . ! . . . . |                                    ",
        "/                           #   | . . . . . . |                                    ",
        "/                           #   | . . . . . . |                                    ",
        "/         - - - -           #   | . . . . . . |                                    ",
        "/         | . . |           # # # # # # # + . . D . . |                            ",
        "/         | < . + # # #       #   | . . . . . . |                                    ",
        "/         - - - -   #         #   | . ? . . . . |                                    ",
        "/                   # # # # # #   - - - - - - - - -                                ",
    ]
    
    # The provided map has inconsistent spacing, this is a corrected version based on the visual structure
    grid = [
        "////////////////////////////////",
        "//|.....|//#########///////////",
        "//|.....|//#       #///////////",
        "//|..g..+### ##### @///////////",
        "//|.....|////#//---+---////////",
        "//-------////#//|.....|////////",
        "/////////////|#//|.!...|////////",
        "/////////////|#//|.....|////////",
        "/////////////|#//|.....|////////",
        "//---//////#//|.....|////////",
        "//|..|//#####+..D..|////////",
        "//|<.+###//#//|.....|////////",
        "//----//#//#//|.?...|////////",
        "///////#####//---------////////",
    ]
    
    # Based on the prompt, let's use a direct interpretation of the provided matrix
    grid = [
        ['/', ' ', ' ', '-', '-', '-', '-', '-', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '],
        ['/', ' ', '|', '.', '.', '.', '.', '|', ' ', ' ', ' ', ' ', ' ', ' ', '#', '#', '#', '#', '#', '#', '#', '#', '#', '#', '#', '#', ' ', ' ', ' ', ' '],
        ['/', ' ', '|', '.', '.', '.', '.', '|', ' ', ' ', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' '],
        ['/', ' ', '|', '.', 'g', '.', '.', '+', '#', '#', '#', '#', '#', '#', '#', '#', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '@', ' ', ' ', ' ', ' '],
        ['/', ' ', '|', '.', '.', '.', '.', '|', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', '-', '-', '-', '+', '-', '-', '-', ' '],
        ['/', ' ', '-', '-', '-', '-', '-', '-', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', '|', '.', '.', '.', '.', '.', '|'],
        ['/', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', '|', '.', '!', '.', '.', '.', '|'],
        ['/', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', '|', '.', '.', '.', '.', '.', '|'],
        ['/', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', '|', '.', '.', '.', '.', '.', '|'],
        ['/', ' ', ' ', ' ', '-', '-', '-', '-', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', '|', '.', '.', '.', '.', '.', '|'],
        ['/', ' ', ' ', ' ', '|', '.', '.', '|', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', '#', '#', '#', '#', '#', '#', '+', '.', '.', 'D', '.', '.', '|'],
        ['/', ' ', ' ', ' ', '|', '<', '.', '+', '#', '#', '#', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', '|', '.', '.', '.', '.', '.', '|'],
        ['/', ' ', ' ', ' ', '-', '-', '-', '-', ' ', ' ', '#', ' ', ' ', ' ', ' ', '#', ' ', ' ', ' ', ' ', ' ', ' ', '|', '.', '?', '.', '.', '.', '|'],
        ['/', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', '#', '#', '#', '#', '#', '#', ' ', ' ', ' ', ' ', ' ', ' ', '-', '-', '-', '-', '-', '-', '-']
    ]


    rows, cols = len(grid), len(grid[0])
    start, goal = None, None
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '@':
                start = (r, c)
            elif grid[r][c] == 'g':
                goal = (r, c)

    walkable_chars = {'.', 'g', '@', '!', '?', '<', '+', '#', ' '}
    
    # Define the dragon's room as a danger zone to avoid
    danger_zone = set()
    for r in range(5, 13):
        for c in range(23, 28):
            danger_zone.add((r,c))

    queue = collections.deque([[start]])
    visited = {start}

    path = None
    while queue:
        current_path = queue.popleft()
        r, c = current_path[-1]

        if (r, c) == goal:
            path = current_path
            break

        for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            nr, nc = r + dr, c + dc

            if 0 <= nr < rows and 0 <= nc < cols and \
               (nr, nc) not in visited and \
               grid[nr][nc] in walkable_chars and \
               (nr, nc) not in danger_zone:
                
                visited.add((nr, nc))
                new_path = list(current_path)
                new_path.append((nr, nc))
                queue.append(new_path)

    if not path:
        print("No safe path found.")
        return

    # Convert path to simplified directions
    directions = []
    last_dir = None
    for i in range(1, len(path)):
        r_prev, c_prev = path[i-1]
        r_curr, c_curr = path[i]
        
        current_dir = None
        if r_curr > r_prev:
            current_dir = 'D'
        elif r_curr < r_prev:
            current_dir = 'U'
        elif c_curr > c_prev:
            current_dir = 'R'
        elif c_curr < c_prev:
            current_dir = 'L'
        
        if current_dir != last_dir:
            directions.append(current_dir)
            last_dir = current_dir
            
    print("".join(directions))

solve_maze()