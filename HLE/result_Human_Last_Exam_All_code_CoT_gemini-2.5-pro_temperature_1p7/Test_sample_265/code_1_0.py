import collections

def solve_bug_climb():
    """
    Solves the climbing bug puzzle by modeling the tower as a graph
    and finding the shortest path using Breadth-First Search (BFS).
    """

    # --- Step 1: Define the properties of the tower and bricks ---
    
    # The repeating color pattern is Red, Red, Blue, Red, Blue, Blue.
    # A brick's color is determined by its index modulo 6.
    # Red bricks correspond to indices 0, 1, and 3 in the pattern.
    RED_MODULOS = {0, 1, 3}

    def is_red(brick_index):
        """Checks if a brick at a given index is red."""
        if brick_index < 0:
            return False
        return (brick_index % 6) in RED_MODULOS

    # --- Step 2: Determine the start and target bricks ---
    start_brick = 0
    circumference = 10.5

    # To be in the same vertical column, the bug must traverse a distance
    # that is a multiple of the circumference and also a whole number.
    # The smallest such positive value is 2 * 10.5 = 21.
    vertical_offset = 2 * circumference
    
    target_brick = -1
    k = 1
    while target_brick == -1:
        potential_target = int(k * vertical_offset)
        if is_red(potential_target):
            target_brick = potential_target
        else:
            k += 1

    # --- Step 3: Define the possible moves (graph edges) ---
    # Moves to adjacent bricks: sideways, and diagonally up/down.
    move_offsets = [-11, -10, -1, 1, 10, 11]

    # --- Step 4: Find the shortest path with Breadth-First Search ---
    # The queue will store tuples of (current_brick, path_taken_list).
    queue = collections.deque([(start_brick, [start_brick])])
    visited = {start_brick}
    
    solution_path = None

    while queue:
        current_brick, path = queue.popleft()

        if current_brick == target_brick:
            solution_path = path
            break

        for offset in move_offsets:
            neighbor_brick = current_brick + offset
            
            # A move is valid if the neighbor is a red brick and hasn't been visited.
            if neighbor_brick not in visited and is_red(neighbor_brick):
                visited.add(neighbor_brick)
                new_path = list(path)
                new_path.append(neighbor_brick)
                queue.append((neighbor_brick, new_path))
    
    # --- Step 5: Output the result ---
    print(f"Start Brick: {start_brick}")
    print(f"Target Brick (first red brick directly above start): {target_brick}")
    
    if solution_path:
        seconds = len(solution_path) - 1
        path_str = " -> ".join(map(str, solution_path))
        print("\nShortest path found:")
        print(f"Path Equation: {path_str}")
        print(f"\nThe minimum number of seconds required is {seconds}.")
    else:
        print("\nAnalysis complete:")
        print("No path exists between the start and target bricks.")
        print("The red bricks form disconnected sections under the given adjacency rules,")
        print("and the bug cannot travel between them.")

# Execute the solver function
solve_bug_climb()