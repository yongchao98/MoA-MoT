import sys

# Increase the recursion limit for deep recursive searches, although not strictly
# necessary for n=10, it's good practice for this type of problem.
sys.setrecursionlimit(2000)

def solve_saw_count(n):
    """
    Calculates the number of n-step self-avoiding walks on a 2D square lattice.

    Args:
        n (int): The number of steps in the walk.

    Returns:
        int: The total number of n-step self-avoiding walks.
    """
    # Possible moves: Up, Down, Left, Right
    moves = [(0, 1), (0, -1), (-1, 0), (1, 0)]
    
    # Path is a set of visited (x, y) tuples for O(1) lookup.
    # The walk starts at the origin (0,0).
    path = {(0, 0)}

    def count_walks(current_pos, steps_left):
        """
        Recursively counts valid self-avoiding walks from the current position.
        """
        # Base case: If we have taken n steps, we found one complete valid walk.
        if steps_left == 0:
            return 1

        total_walks = 0
        x, y = current_pos

        # Explore all four possible directions.
        for dx, dy in moves:
            next_pos = (x + dx, y + dy)
            
            # If the next position has not been visited yet...
            if next_pos not in path:
                # ...take the step: add the new position to the path.
                path.add(next_pos)
                
                # Recursively continue the walk from the new position.
                total_walks += count_walks(next_pos, steps_left - 1)
                
                # Backtrack: remove the position to explore other branches.
                path.remove(next_pos)
        
        return total_walks

    # Start the calculation from the origin (0,0) with n steps to take.
    return count_walks((0, 0), n)

if __name__ == "__main__":
    n = 10
    result = solve_saw_count(n)
    # The problem asks to output each number in the final equation.
    # The equation is a(10) = result. The numbers are 10 and the result.
    print(f"a({n}) = {result}")
