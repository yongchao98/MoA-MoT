import sys

# It's good practice to increase the recursion limit for deep recursion,
# though it may not be strictly necessary for n=10.
sys.setrecursionlimit(2000)

def count_walks(steps_left, current_pos, visited):
    """
    Recursively counts the number of self-avoiding walks.

    :param steps_left: The number of steps remaining to take.
    :param current_pos: The current (x, y) coordinates on the lattice.
    :param visited: A set of (x, y) tuples representing visited points.
    :return: The total number of valid walks from the current state.
    """
    # Base case: If there are no steps left, we have successfully found one valid walk.
    if steps_left == 0:
        return 1

    total_walks = 0
    x, y = current_pos
    # Possible moves: North, South, East, West
    directions = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    # Explore each possible direction
    for dx, dy in directions:
        next_pos = (x + dx, y + dy)
        # Check if the next position has already been visited
        if next_pos not in visited:
            # If not visited, add it to the visited set
            visited.add(next_pos)
            # Recursively call the function for the next step
            total_walks += count_walks(steps_left - 1, next_pos, visited)
            # Backtrack: remove the position to explore other paths
            visited.remove(next_pos)

    return total_walks

def solve_saw_count(n):
    """
    Initializes and solves the self-avoiding walk problem for n steps.
    """
    # The walk starts at the origin (0,0)
    start_pos = (0, 0)
    # The initial set of visited points contains only the origin
    visited = {start_pos}
    
    # Calculate the total number of n-step self-avoiding walks
    result = count_walks(n, start_pos, visited)
    
    # Print the final equation as requested
    print(f"a({n}) = {result}")

if __name__ == '__main__':
    # The problem asks for a(10)
    number_of_steps = 10
    solve_saw_count(number_of_steps)
