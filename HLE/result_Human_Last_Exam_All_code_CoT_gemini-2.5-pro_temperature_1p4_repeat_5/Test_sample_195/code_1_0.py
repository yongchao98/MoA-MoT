import heapq

def solve_dungeon_path():
    """
    Finds the least dangerous path from the adventurer '@' to the gold 'g'
    by using the A* pathfinding algorithm, avoiding the Red Dragon 'D'.
    """
    # The map layout represented as a grid.
    grid = [
        "                                         ",
        "   |.......|        ######################   ",
        "   |.......|        #                    #   ",
        "   |.g..+##################@           #   ",
        "   |.......|        #       ---+---        #   ",
        "   ---------        #       |.......|        ",
        "                    #       |..!....|        ",
        "                    #       |.......|        ",
        "                    #       |.......|        ",
        "       ----         #       |.......|        ",
        "       |..|     #######+|...D...|        ",
        "       |<.+#########  #       |.......|        ",
        "       ----         #       |...?...|        ",
        "                    ####### |.......|        ",
        "                            ---------        "
    ]

    walls = {'#', '|', '-'}
    start, end, dragon = None, None, None

    # Parse the grid to find key locations
    for r, row in enumerate(grid):
        for c, char in enumerate(row):
            if char == '@':
                start = (r, c)
            elif char == 'g':
                end = (r, c)
            elif char == 'D':
                dragon = (r, c)

    def get_cost(pos):
        """
        Determines the cost of moving to a tile.
        The danger zone is the large chamber containing the dragon.
        """
        r, c = pos
        # Define the bounding box for the dragon's chamber
        if 5 <= r <= 14 and 24 <= c <= 34:
            return 100  # High cost for danger zone
        return 1 # Normal cost for other tiles

    def heuristic(a, b):
        """Calculates Manhattan distance as the heuristic for A*."""
        return abs(a[0] - b[0]) + abs(a[1] - b[1])

    def find_path():
        """A* pathfinding algorithm implementation."""
        open_set = []
        # Heap items are: (f_score, g_score, position, path_so_far)
        heapq.heappush(open_set, (0 + heuristic(start, end), 0, start, []))
        closed_set = set()

        while open_set:
            f_score, g_score, current, path = heapq.heappop(open_set)

            if current == end:
                return path + [current]

            if current in closed_set:
                continue
            closed_set.add(current)
            
            r, c = current
            # Explore neighbors (Up, Down, Left, Right)
            for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                neighbor = (r + dr, c + dc)
                nr, nc = neighbor

                # Check bounds and walls
                if not (0 <= nr < len(grid) and 0 <= nc < len(grid[0]) and grid[nr][nc] not in walls):
                    continue
                
                if neighbor in closed_set:
                    continue

                new_g_score = g_score + get_cost(neighbor)
                new_f_score = new_g_score + heuristic(neighbor, end)
                
                # Check if this neighbor is already in the open set with a lower f_score
                in_open_set = False
                for item in open_set:
                    if item[2] == neighbor and item[0] <= new_f_score:
                        in_open_set = True
                        break
                
                if not in_open_set:
                    heapq.heappush(open_set, (new_f_score, new_g_score, neighbor, path + [current]))
        return None # No path found

    path_coords = find_path()
    
    if not path_coords:
        print("No path could be found.")
        return

    # Convert coordinate path to a sequence of moves
    moves = []
    for i in range(1, len(path_coords)):
        prev_r, prev_c = path_coords[i-1]
        curr_r, curr_c = path_coords[i]
        if curr_r > prev_r: moves.append('D')
        elif curr_r < prev_r: moves.append('U')
        elif curr_c > prev_c: moves.append('R')
        elif curr_c < prev_c: moves.append('L')

    # Compress the move sequence as per instructions
    if not moves:
        print("Already at the destination.")
        return
        
    compressed_moves = [moves[0]]
    for move in moves[1:]:
        if move != compressed_moves[-1]:
            compressed_moves.append(move)

    final_path_str = "".join(compressed_moves)
    print("The least dangerous path is:")
    print(final_path_str)
    
    # Printing the final answer in the requested format
    print("\n<<<" + final_path_str + ">>>")

solve_dungeon_path()