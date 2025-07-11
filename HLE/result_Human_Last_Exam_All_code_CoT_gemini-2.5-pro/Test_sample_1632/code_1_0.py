def count_saw(x, y, n, visited):
    """
    Recursively counts the number of self-avoiding walks.

    Args:
        x: The current x-coordinate.
        y: The current y-coordinate.
        n: The number of steps remaining in the walk.
        visited: A set of (x, y) tuples representing the path taken so far.

    Returns:
        The number of n-step self-avoiding walks from the current position.
    """
    # Base case: If no steps are left, we have found one complete walk.
    if n == 0:
        return 1

    count = 0
    # Define the four possible moves: North, South, East, West.
    moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    # Explore each possible move.
    for dx, dy in moves:
        nx, ny = x + dx, y + dy
        
        # If the next position has not been visited yet.
        if (nx, ny) not in visited:
            # Add the new position to the visited set for the recursive call.
            visited.add((nx, ny))
            # Recursively call the function for the next step.
            count += count_saw(nx, ny, n - 1, visited)
            # Backtrack: remove the new position to explore other paths from the current state.
            visited.remove((nx, ny))
    
    return count

def main():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a Manhattan lattice.
    """
    N = 10
    
    # We use symmetry to speed up the calculation. The total number of walks is
    # 4 * (number of walks starting with a step in a specific direction, e.g., North).
    # We calculate the number of walks starting with a step from (0,0) to (0,1).
    # This means the walk has already taken 1 step.
    
    # The path starts at (0,0) and takes the first step to (0,1).
    initial_visited = {(0, 0), (0, 1)}
    
    # We need to find walks of length N-1 starting from (0,1).
    steps_remaining = N - 1
    start_x, start_y = 0, 1
    
    # Calculate the number of walks starting with a step North.
    walks_starting_north = count_saw(start_x, start_y, steps_remaining, initial_visited)
    
    # The total number of walks is 4 times this value due to symmetry.
    total_walks = 4 * walks_starting_north
    
    # Print the final equation with all its components.
    print(f"a(10) = 4 * {walks_starting_north} = {total_walks}")

if __name__ == "__main__":
    main()