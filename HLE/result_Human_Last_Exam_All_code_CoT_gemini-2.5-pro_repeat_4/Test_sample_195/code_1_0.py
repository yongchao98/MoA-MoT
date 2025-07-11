import collections

def solve_maze():
    """
    This function finds the least dangerous path from the adventurer '@' to the gold 'g'.

    My plan is as follows:
    1.  **Represent the Map**: I will store the provided ASCII map as a list of strings.
    2.  **Analyze the Rules**:
        - Walls are '/', '-', '|'. These are impassable.
        - The goal is 'g' (gold). The start is '@' (adventurer).
        - A '#' represents an "Unlit hallway", which is traversable.
        - The main danger is the Red Dragon 'D'. The "least dangerous" path is one that avoids the dragon.
    3.  **Find the Path**: I will use a Breadth-First Search (BFS) algorithm. BFS is ideal for finding the shortest path in an unweighted grid like this. The shortest path will also be the least dangerous because it avoids the dragon's area entirely, which is on a completely different route.
    4.  **Determine Directions**: Once the BFS finds the path (as a sequence of coordinates), I will convert this into a series of moves (U, D, L, R).
    5.  **Simplify and Print**: As per the instructions, I will simplify the sequence of moves (e.g., 'LLL' becomes 'L') and print the final result.
    """
    # The map is represented as a list of strings. I've corrected the alignment
    # from the prompt's visual layout to create a usable grid.
    MAP = [
        "           - - - - - - -                    ",
        "           | . . . . |                      ##########",
        "           | . . . . |                      #          #",
        "           | . g . . + #################### #        @ #",
        "           - - - - - -          #          #   - - - + - - - ",
        "                              #          #   | . . . . . |",
        "                              #          #   | . ! . . . |",
        "                              #          #   | . . . . . |",
        "                              #          #   | . . . . . |",
        "                 - - - -      #          #   | . . . . . |",
        "                 | . . |      #      #####+ . . D . . |",
        "                 | < . + ###  #      #      | . . . . . |",
        "                 - - - -   #  #      #      | . ? . . . |",
        "                           #####      #      - - - - - - - - "
    ]

    # Find the starting and ending coordinates
    try:
        start_pos = next((r, c) for r, row in enumerate(MAP) for c, char in enumerate(row) if char == '@')
        goal_pos = next((r, c) for r, row in enumerate(MAP) for c, char in enumerate(row) if char == 'g')
    except StopIteration:
        print("Error: Start ('@') or goal ('g') not found in the map.")
        return

    # Define obstacles based on the key
    obstacles = {'/', '-', '|'}

    # Initialize BFS queue with (position, path_list)
    queue = collections.deque([(start_pos, [start_pos])])
    visited = {start_pos}
    path_to_goal = None

    while queue:
        (r, c), path = queue.popleft()

        if (r, c) == goal_pos:
            path_to_goal = path
            break

        # Explore neighbors (Up, Down, Left, Right)
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nr, nc = r + dr, c + dc

            # Check if neighbor is within bounds, not an obstacle, and not visited
            if (0 <= nr < len(MAP) and 0 <= nc < len(MAP[nr]) and
                    MAP[nr][nc] not in obstacles and (nr, nc) not in visited):
                visited.add((nr, nc))
                new_path = path + [(nr, nc)]
                queue.append(((nr, nc), new_path))
    
    if not path_to_goal:
        print("No path could be found.")
        return

    # Convert the list of coordinates into a list of direction characters
    directions = []
    for i in range(1, len(path_to_goal)):
        prev_r, prev_c = path_to_goal[i-1]
        curr_r, curr_c = path_to_goal[i]
        if curr_r > prev_r: directions.append('D')
        elif curr_r < prev_r: directions.append('U')
        elif curr_c > prev_c: directions.append('R')
        elif curr_c < prev_c: directions.append('L')

    # Simplify the list of directions (e.g., L, L, L -> L)
    if not directions:
        # This case handles if start and goal are the same
        return
        
    simplified_directions = [directions[0]]
    for direction in directions[1:]:
        if direction != simplified_directions[-1]:
            simplified_directions.append(direction)
    
    # Print the final simplified path
    print("".join(simplified_directions))

solve_maze()