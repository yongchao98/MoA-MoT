import sys

def count_self_avoiding_walks(steps_left, current_pos, visited):
    """
    Recursively counts the number of self-avoiding walks.

    Args:
        steps_left (int): The number of steps remaining to take.
        current_pos (tuple): The current (x, y) coordinate of the walk.
        visited (set): A set of (x, y) tuples representing the visited points
                       in the current path.

    Returns:
        int: The number of valid self-avoiding walks that can be formed from
             the current state.
    """
    # Base case: If there are no steps left, we have found one valid walk.
    if steps_left == 0:
        return 1

    # Get the current coordinates.
    x, y = current_pos
    
    # Initialize a counter for the number of walks from this point.
    walk_count = 0
    
    # Define the four possible moves on a Manhattan lattice.
    moves = [(0, 1), (0, -1), (1, 0), (-1, 0)] # Up, Down, Right, Left

    # Explore each possible move.
    for dx, dy in moves:
        next_pos = (x + dx, y + dy)
        
        # Check if the next position has already been visited.
        if next_pos not in visited:
            # If not visited, explore this new path.
            # 1. Add the new position to the visited set.
            visited.add(next_pos)
            
            # 2. Recursively call the function for the next step.
            walk_count += count_self_avoiding_walks(steps_left - 1, next_pos, visited)
            
            # 3. Backtrack: Remove the new position to explore other paths.
            visited.remove(next_pos)
            
    return walk_count

def main():
    """
    Calculates and prints the number of 10-step self-avoiding walks.
    """
    # The number of steps for the walk, n.
    n = 10
    
    # A self-avoiding walk starts at an origin point, e.g., (0, 0).
    start_position = (0, 0)
    
    # The set of visited points initially contains only the starting position.
    # Using a set provides fast O(1) average time complexity for lookups.
    initial_visited = {start_position}
    
    # Calculate the total number of walks by starting the recursion.
    total_walks = count_self_avoiding_walks(n, start_position, initial_visited)
    
    print("Let a(n) be the number of n-step self-avoiding walks on a Manhattan lattice.")
    print(f"For n = {n}, the number of walks is:")
    print(f"a({n}) = {total_walks}")

if __name__ == "__main__":
    # For deep recursion, it might be necessary to increase the recursion limit,
    # but n=10 is well within the default limit for Python.
    # sys.setrecursionlimit(20000)
    main()