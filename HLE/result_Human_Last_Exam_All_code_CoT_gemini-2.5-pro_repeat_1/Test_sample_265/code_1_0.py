import math
from collections import deque

def solve_bug_path():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    # The color pattern is R, R, B, R, B, B, repeating every 6 bricks.
    # Red bricks are at relative positions 0, 1, and 3 in the pattern.
    def is_red(brick_index):
        """Checks if a brick is red based on its index in the coil."""
        return brick_index % 6 in [0, 1, 3]

    # A brick at index `i` is adjacent to its neighbors along the coil (i-1, i+1)
    # and the bricks it touches in the rows above and below. With a circumference
    # of 10.5, these are the integer bricks surrounding i +/- 10.5.
    def get_neighbors(brick_index):
        """Returns the indices of all physically adjacent bricks."""
        c = 10.5  # Circumference
        return [
            brick_index - 1,
            brick_index + 1,
            math.floor(brick_index - c),
            math.ceil(brick_index - c),
            math.floor(brick_index + c),
            math.ceil(brick_index + c),
        ]

    # A goal brick must be above the start (index > 0) and in the same
    # vertical column (index is a multiple of 21).
    def is_goal(brick_index):
        """Checks if a brick is a valid goal destination."""
        if brick_index > 0 and brick_index % 21 == 0:
            return True
        return False

    # Initialize the BFS
    start_node = 0
    # The queue stores tuples of (current_node, path_taken_to_get_here).
    queue = deque([(start_node, [start_node])])
    # 'visited' stores nodes already queued to avoid cycles and redundant work.
    visited = {start_node}

    print("Searching for the shortest path...")
    
    while queue:
        current_node, path = queue.popleft()

        # If the current node is a valid goal, we've found the shortest path.
        if is_goal(current_node) and is_red(current_node):
            seconds = len(path) - 1
            print("\n--- Solution Found ---")
            print(f"The bug can reach a brick at index {current_node}, which is vertically above the start.")
            print(f"The minimum time required is {seconds} seconds.")
            print("\nThe shortest path is:")
            path_str = " -> ".join(map(str, path))
            print(f"{path_str} = {seconds}")
            return seconds

        # Explore the neighbors of the current brick
        for neighbor in get_neighbors(current_node):
            if neighbor not in visited and is_red(neighbor):
                visited.add(neighbor)
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))
    
    print("No path found.") # This should not be reached for this problem
    return None

# Execute the search and print the result.
final_answer = solve_bug_path()
if final_answer is not None:
    print(f"\n<<< {final_answer} >>>")
