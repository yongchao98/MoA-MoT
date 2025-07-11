def solve_self_avoiding_walk():
    """
    Calculates a(n), the number of n-step self-avoiding walks on a Manhattan lattice.
    This solution uses a recursive backtracking algorithm with memoization to count the walks efficiently.
    """
    TARGET_LENGTH = 10
    memo = {}

    def count_walks(length, pos, visited):
        """
        Recursively counts valid self-avoiding walks.

        Args:
            length: The current length of the walk.
            pos: The current (x, y) coordinates of the head of the walk.
            visited: A set of (x, y) tuples of all points in the walk.

        Returns:
            The number of valid self-avoiding walks of TARGET_LENGTH starting from the given state.
        """
        # Base case: if the walk has reached the target length, we've found one valid path.
        if length == TARGET_LENGTH:
            return 1

        # To use memoization, we create a canonical key representing the path's shape
        # relative to its current endpoint. This makes the cache useful for paths that are
        # translations of each other.
        px, py = pos
        relative_visited = frozenset((vx - px, vy - py) for vx, vy in visited)
        
        memo_key = (length, relative_visited)
        if memo_key in memo:
            return memo[memo_key]

        count = 0
        x, y = pos
        
        # Explore the four possible directions on the Manhattan lattice.
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            next_pos = (x + dx, y + dy)
            
            # If the next position has not been visited, continue the walk from there.
            if next_pos not in visited:
                # Create a new visited set for the recursive call.
                new_visited = visited | {next_pos}
                count += count_walks(length + 1, next_pos, new_visited)

        # Store the result in the cache before returning.
        memo[memo_key] = count
        return count

    # Start the walk at the origin (0,0). The initial state is a walk of length 0.
    initial_pos = (0, 0)
    initial_visited = {initial_pos}
    result = count_walks(0, initial_pos, initial_visited)

    print(f"a({TARGET_LENGTH}) = {result}")

if __name__ == '__main__':
    solve_self_avoiding_walk()