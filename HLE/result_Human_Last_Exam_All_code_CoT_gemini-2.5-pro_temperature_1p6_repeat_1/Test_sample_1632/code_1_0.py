import sys

# It's good practice to increase the recursion limit for deep recursion,
# though n=10 is shallow enough for default limits.
sys.setrecursionlimit(2000)

def count_walks(current_pos, steps_left, visited):
    """
    Recursively counts the number of self-avoiding walks from a given position.

    Args:
        current_pos (tuple): The (x, y) coordinates of the current position.
        steps_left (int): The number of steps remaining in the walk.
        visited (set): A set of (x, y) tuples representing visited positions.

    Returns:
        int: The number of valid self-avoiding walks from the current state.
    """
    # Base case: If there are no steps left, we have successfully found one complete walk.
    if steps_left == 0:
        return 1

    x, y = current_pos
    total_walks = 0

    # The four possible moves on a Manhattan lattice.
    moves = [(0, 1),  # Up
             (0, -1), # Down
             (1, 0),  # Right
             (-1, 0)] # Left

    # Explore each possible move.
    for dx, dy in moves:
        next_pos = (x + dx, y + dy)
        # A move is valid if the next position has not been visited yet.
        if next_pos not in visited:
            # Mark the new position as visited for the subsequent recursive calls.
            visited.add(next_pos)
            # Recurse with one less step.
            total_walks += count_walks(next_pos, steps_left - 1, visited)
            # Backtrack: Un-mark the position so it can be used by other paths.
            visited.remove(next_pos)

    return total_walks

def main():
    """
    Calculates a(n), the number of n-step self-avoiding walks.
    """
    n = 10

    if n == 0:
        result = 1
        print(f"a(0) = {result}")
        return

    # Optimization: Exploit the symmetry of the grid.
    # The number of walks starting with a move in any of the 4 directions is the same.
    # We calculate the number of walks starting with a move to the right,
    # and then multiply by 4.

    # The walk starts at (0,0). The first step is fixed to (1,0).
    first_step_pos = (1, 0)
    
    # The set of visited points must include the origin and the first step.
    visited_path = {(0, 0), first_step_pos}

    # We now need to find the number of (n-1)-step walks from the new position.
    steps_remaining = n - 1
    
    walks_per_branch = count_walks(first_step_pos, steps_remaining, visited_path)
    
    # The total number of walks is 4 times the walks in one branch.
    total_walks = 4 * walks_per_branch

    print(f"Let a(n) be the number of {n}-step self-avoiding walks.")
    print(f"We can find a({n}) by calculating the walks for one initial direction and multiplying by 4 due to symmetry.")
    print(f"a({n}) = 4 * {walks_per_branch} = {total_walks}")

if __name__ == "__main__":
    main()