import collections

def solve_tower_puzzle():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).

    This function models the tower as a graph where bricks are nodes and adjacency
    is defined by the coil structure. It finds the shortest path from the start
    brick (0) to a target brick (in the same vertical column) using only red bricks.
    """

    # A helper function to determine if a brick at a given index 'n' is red.
    # The pattern is RRBRBB, repeating every 6 bricks.
    # Red bricks correspond to indices 0, 1, and 3 within the pattern.
    def is_red(n):
        return n % 6 in {0, 1, 3}

    # Initialize the BFS queue with the starting state.
    # Each item in the queue is a tuple: (current_brick_index, path_taken_so_far).
    queue = collections.deque([(0, [0])])

    # A set to keep track of visited bricks to prevent cycles and redundant computations.
    visited = {0}

    # The bug can move to 6 adjacent bricks based on the coil structure.
    # Moves of +/-1 are along the coil.
    # Moves of +/-10 and +/-11 are to the staggered bricks above/below.
    moves = [-11, -10, -1, 1, 10, 11]

    # Main BFS loop
    while queue:
        current_brick, path = queue.popleft()

        # Check if the current brick is a target.
        # A target is in the same vertical column, so its index must be a
        # positive multiple of 21 (LCM of 1 and 10.5).
        if current_brick > 0 and current_brick % 21 == 0:
            # We found the shortest path to a target.
            # The number of seconds is the number of steps in the path.
            num_seconds = len(path) - 1
            equation_str = " -> ".join(map(str, path))

            print("The final equation representing the shortest path is:")
            print(f"{equation_str} = {num_seconds} seconds")
            return

    # This part of the code should not be reached if a solution exists.
    print("No solution found.")


if __name__ == '__main__':
    solve_tower_puzzle()