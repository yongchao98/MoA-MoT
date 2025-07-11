import sys

def solve_saw_count():
    """
    Calculates a(n), the number of n-step self-avoiding walks on a Manhattan lattice,
    using a recursive backtracking algorithm with a symmetry optimization.
    """
    n = 10

    # The recursive function to count walks from a given state.
    def count_walks_recursive(steps_left, current_pos, visited_path):
        """
        Args:
            steps_left (int): The number of steps remaining to take.
            current_pos (tuple): The current (x, y) coordinates.
            visited_path (set): A set of visited (x, y) coordinates.
        """
        # Base case: If no steps are left, we have found one valid walk.
        if steps_left == 0:
            return 1

        x, y = current_pos
        total = 0
        
        # Explore the four possible moves (Up, Down, Right, Left).
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            next_pos = (x + dx, y + dy)

            # Check if the next position has already been visited.
            if next_pos not in visited_path:
                # If not visited, take the step.
                visited_path.add(next_pos)
                
                # Recursively call the function for the next step.
                total += count_walks_recursive(steps_left - 1, next_pos, visited_path)
                
                # Backtrack: undo the move for exploring other branches.
                visited_path.remove(next_pos)
        
        return total

    # --- Calculation using Symmetry ---
    # We calculate walks for n-1 steps, starting from a point adjacent to the origin,
    # and then multiply by 4.
    
    # We want a(10), so we take 1 initial step and then need to find paths of length 9.
    steps_after_first = n - 1
    
    # Assume the first step is from (0,0) to (1,0).
    first_step_pos = (1, 0)
    
    # The initial path already contains the origin and the first step.
    initial_path = {(0, 0), first_step_pos}

    # Calculate the number of valid walks for this one starting direction.
    walks_in_one_direction = count_walks_recursive(steps_after_first, first_step_pos, initial_path)

    # The total number of walks is 4 times this value due to symmetry.
    total_walks = 4 * walks_in_one_direction
    
    # --- Output ---
    print(f"To find a({n}), the number of {n}-step self-avoiding walks:")
    print("We calculate the number of walks that start with a step in one direction and multiply by 4.")
    print(f"Number of {n-1}-step walks starting from (1,0) (given a path from (0,0)): {walks_in_one_direction}")
    print(f"Final Equation: a({n}) = 4 * {walks_in_one_direction}")
    print(f"Result: a({n}) = {total_walks}")
    
    # The problem requires the final answer to be wrapped in <<<>>>
    # The execution of this function prints the requested information.
    # The final answer will be provided after the code block.

solve_saw_count()