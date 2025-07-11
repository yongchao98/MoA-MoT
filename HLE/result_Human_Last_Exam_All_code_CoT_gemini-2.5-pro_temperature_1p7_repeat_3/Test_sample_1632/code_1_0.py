import sys

def solve_saw_count(n):
    """
    Calculates a(n), the number of n-step self-avoiding walks on a Manhattan lattice.
    
    This function uses a recursive backtracking algorithm to explore all possible paths.
    
    Args:
        n (int): The number of steps in the walk.
        
    Returns:
        int: The total number of n-step self-avoiding walks.
    """
    # Increase recursion limit if needed for larger n, though not necessary for n=10.
    # sys.setrecursionlimit(max(sys.getrecursionlimit(), n + 5))

    # A set to store visited coordinates for the current path for O(1) lookups.
    visited = set()
    # The walk starts at the origin.
    visited.add((0, 0))

    def count_walks_recursive(x, y, steps_left):
        """
        Recursively counts valid self-avoiding walks from the current position.
        
        Args:
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            steps_left (int): The number of steps remaining in the walk.
            
        Returns:
            int: The number of valid walks from this state.
        """
        # Base case: If no steps are left, we have found one complete walk.
        if steps_left == 0:
            return 1

        total_walks = 0
        # The four possible moves on a Manhattan lattice: up, down, left, right.
        moves = [(0, 1), (0, -1), (1, 0), (-1, 0)]

        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy
            next_point = (next_x, next_y)

            # Check if the next point would violate the self-avoiding condition.
            if next_point not in visited:
                # Explore from the new point.
                visited.add(next_point)
                total_walks += count_walks_recursive(next_x, next_y, steps_left - 1)
                # Backtrack: remove the point to allow other paths to use it.
                visited.remove(next_point)
        
        return total_walks

    # Start the walk from the origin (0, 0) with n steps.
    result = count_walks_recursive(0, 0, n)
    return result

if __name__ == '__main__':
    # The number of steps for the walk.
    num_steps = 10
    
    # Calculate the result.
    a_10 = solve_saw_count(num_steps)
    
    # The problem asks for the equation, which is a(n) = result.
    # We print each number involved in this final equation.
    print(f"Let a(n) be the number of n-step self-avoiding walks.")
    print(f"For n = {num_steps}:")
    print(f"a({num_steps}) = {a_10}")
