import sys

def solve_saw(n):
    """
    Calculates the number a(n) of n-step self-avoiding walks on a square lattice.
    
    The method uses a backtracking algorithm. To optimize, it leverages the
    symmetry of the lattice: the total number of walks is 4 times the number
    of walks that start with a specific first step (e.g., to the East).
    """
    if n == 0:
        return 1

    # Start at (0,0). The first step is fixed to (1,0) to use symmetry.
    # The 'visited' set starts with the origin and the first step's destination.
    visited = {(0, 0), (1, 0)}

    def count_recursive(x, y, steps_left):
        """
        Recursively counts the number of valid walks from the current position (x, y).
        'visited' is a closure variable from the parent function's scope.
        """
        # Base case: if no steps are left, we have found one complete valid walk.
        if steps_left == 0:
            return 1

        total = 0
        
        # Define the four possible moves from the current position (x,y).
        moves = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

        for next_x, next_y in moves:
            next_pos = (next_x, next_y)
            if next_pos not in visited:
                # If the next position is valid (not visited), we explore it.
                visited.add(next_pos)
                
                # Recursively count the walks from this new position.
                total += count_recursive(next_x, next_y, steps_left - 1)
                
                # Backtrack: un-mark the position to allow other paths to use it.
                visited.remove(next_pos)
        
        return total

    # The recursion starts from (1,0) with n-1 steps remaining.
    count_from_one_direction = count_recursive(1, 0, n - 1)
    
    # The total number of walks is 4 times the count from one direction due to symmetry.
    return 4 * count_from_one_direction

# Set the number of steps for the walk.
n_steps = 10
# Calculate the result using the function.
result = solve_saw(n_steps)

# Print the final result in the requested format "a(n) = result".
print(f"a({n_steps}) = {result}")
