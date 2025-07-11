import sys

def solve_a_n():
    """
    Calculates and prints a(n), the number of n-step self-avoiding walks 
    on a 2D square lattice (Manhattan lattice).
    """
    # The number of steps for the walk.
    n = 10
    
    # We use a set to keep track of visited coordinates for efficient lookups.
    # The walk starts at the origin (0, 0).
    visited = {(0, 0)}
    
    # Start the recursive counting process from the origin.
    result = _count_recursive(n, (0, 0), visited)
    
    # As requested, the script prints the final result. The equation format
    # "a(n) = result" includes both numbers from the problem statement and the solution.
    print(f"a({n}) = {result}")

def _count_recursive(steps_left, current_pos, visited):
    """
    A recursive helper function that explores paths and counts the valid ones.
    
    :param steps_left: The number of steps remaining in the walk.
    :param current_pos: The current (x, y) coordinate in the walk.
    :param visited: A set of (x, y) coordinates already visited in the current path.
    :return: The number of valid self-avoiding walks from the current state.
    """
    
    # Base case: If we have completed all steps, we have found one valid walk.
    if steps_left == 0:
        return 1

    count = 0
    x, y = current_pos
    
    # Define the four possible moves on a Manhattan lattice (North, East, South, West).
    moves = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

    # Iterate through each potential next move.
    for next_pos in moves:
        # A move is valid only if the destination has not been visited before.
        if next_pos not in visited:
            # If the move is valid:
            # 1. Mark the new position as visited.
            visited.add(next_pos)
            
            # 2. Recursively call the function for the next step and add the
            #    number of walks found to the total count.
            count += _count_recursive(steps_left - 1, next_pos, visited)
            
            # 3. Backtrack: remove the position from the visited set. This is
            #    essential for the algorithm to explore all distinct paths correctly.
            visited.remove(next_pos)
            
    return count

if __name__ == "__main__":
    # To prevent potential recursion depth issues with larger n, though not expected for n=10.
    sys.setrecursionlimit(2000)
    solve_a_n()