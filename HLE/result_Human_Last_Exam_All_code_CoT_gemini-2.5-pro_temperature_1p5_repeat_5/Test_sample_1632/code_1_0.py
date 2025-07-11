import sys

def solve_saw_count():
    """
    Calculates a(10), the number of 10-step self-avoiding walks on a 
    Manhattan (2D square) lattice.

    The method uses a recursive backtracking algorithm. A walk starts at (0,0) 
    and explores all possible paths of a given length, ensuring no point is 
    visited more than once.

    To optimize the calculation, we use the symmetry of the square lattice. The 
    total number of walks is 4 times the number of walks that start with a 
    step in a specific direction (e.g., from (0,0) to (1,0)). The recursion 
    then calculates the number of remaining 9-step walks starting from (1,0).
    """
    # The desired number of steps for the walk.
    n = 10
    
    # We increase the recursion limit as a safeguard for deep recursive calls,
    # though Python's default limit is usually sufficient for n=10.
    sys.setrecursionlimit(2000)

    def count_walks_recursive(pos, visited, steps_left):
        """
        Recursively counts the number of self-avoiding walks from the current state.

        Args:
            pos (tuple): The current (x, y) position.
            visited (set): A set of visited (x, y) position tuples.
            steps_left (int): The number of steps remaining to take.
        """
        # Base case: if no steps are left, we have found one complete walk.
        if steps_left == 0:
            return 1

        x, y = pos
        count = 0
        
        # Explore the 4 possible neighbor positions.
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_pos = (x + dx, y + dy)
            
            # If the neighbor has not been visited, explore from there.
            if next_pos not in visited:
                visited.add(next_pos)
                count += count_walks_recursive(next_pos, visited, steps_left - 1)
                visited.remove(next_pos)  # Backtrack to explore other paths.
        
        return count

    # Using the symmetry optimization:
    # We calculate the number of walks for one initial direction and multiply by 4.
    # Initial state for recursion:
    # - Start at (0,0), take the first step to (1,0).
    # - The path now contains (0,0) and (1,0).
    # - There are n-1 = 9 steps left to take.
    
    start_pos = (0, 0)
    first_step_pos = (1, 0)
    initial_visited = {start_pos, first_step_pos}
    steps_remaining = n - 1
    
    # Calculate the number of walks for one initial direction.
    walks_one_direction = count_walks_recursive(
        first_step_pos, initial_visited, steps_remaining
    )

    # The total number of walks is 4 times this amount.
    total_walks = 4 * walks_one_direction
    
    # Print the final result in an equation format.
    print(f"To find a(10), we count the walks starting in one direction and multiply by 4.")
    print(f"Number of 10-step walks starting with a step to the right: {walks_one_direction}")
    print(f"a(10) = 4 * {walks_one_direction} = {total_walks}")

solve_saw_count()