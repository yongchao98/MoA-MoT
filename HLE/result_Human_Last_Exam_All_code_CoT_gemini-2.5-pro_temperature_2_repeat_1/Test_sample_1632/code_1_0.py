def count_walks_recursive(x, y, steps_left, visited):
    """
    A recursive function to count the number of self-avoiding walks.
    - x, y: current coordinates of the walk.
    - steps_left: number of steps remaining in the walk.
    - visited: a set of (x,y) tuples representing the path so far.
    """
    # Base case: if there are no steps left, we have found one valid walk.
    if steps_left == 0:
        return 1

    count = 0
    # Explore neighbors: up, down, left, right
    moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    for dx, dy in moves:
        next_x, next_y = x + dx, y + dy

        # If the neighbor has not been visited, we can move there.
        if (next_x, next_y) not in visited:
            # Add the new point to the visited set for the recursive call.
            visited.add((next_x, next_y))
            # Recursively count walks from the new point.
            count += count_walks_recursive(next_x, next_y, steps_left - 1, visited)
            # Backtrack: remove the point to explore other paths.
            visited.remove((next_x, next_y))
            
    return count

def main():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a 
    Manhattan lattice.
    """
    n = 10
    
    print(f"Let a(n) be the number of n-step self-avoiding walks on a Manhattan lattice.")
    print(f"To find a({n}), we assume the walk starts at (0,0).")
    print(f"By symmetry, a({n}) is 4 times the number of walks that start by moving in one specific direction (e.g., to the right).")

    # A walk of n steps starting at (0,0) and moving right first
    # is equivalent to an (n-1) step walk from (1,0) that avoids (0,0).
    steps_remaining = n - 1
    
    # The path already contains the starting point (0,0) and the first step to (1,0).
    start_point = (0, 0)
    first_step_pos = (1, 0)
    initial_visited = {start_point, first_step_pos}

    # Calculate the number of walks for one initial direction.
    walks_one_direction = count_walks_recursive(
        first_step_pos[0], 
        first_step_pos[1], 
        steps_remaining, 
        initial_visited
    )

    # The total number of walks is 4 times this value.
    total_walks = 4 * walks_one_direction

    print(f"\nThe number of {n}-step walks starting with a move to the right is: {walks_one_direction}")
    print(f"The final calculation is:")
    print(f"a({n}) = 4 * {walks_one_direction} = {total_walks}")

if __name__ == '__main__':
    main()