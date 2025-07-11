import sys

def solve():
    """
    Calculates a(n), the number of n-step self-avoiding walks on a Manhattan lattice.
    This script specifically solves for n=10.
    """
    # The desired number of steps in the walk.
    n = 10
    
    # Increase recursion limit if needed, though default is usually fine for n=10.
    # sys.setrecursionlimit(max(sys.getrecursionlimit(), n + 5))

    # Using a mutable list to hold the counter so it can be modified
    # within the recursive function's scope.
    count = [0]

    def find_walks(path, x, y, length):
        """
        Recursively finds and counts self-avoiding walks.
        
        Args:
            path (set): A set of (x, y) tuples representing visited points.
            x (int): The current x-coordinate.
            y (int): The current y-coordinate.
            length (int): The current length of the walk (number of steps taken).
        """
        # Base case: If the walk has reached the desired length, count it.
        if length == n:
            count[0] += 1
            return

        # The four possible moves on a Manhattan lattice.
        moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]

        # Explore each possible move.
        for dx, dy in moves:
            next_x, next_y = x + dx, y + dy

            # Check if the next point is already in the path.
            if (next_x, next_y) not in path:
                # If not, add it to the path and recurse.
                path.add((next_x, next_y))
                find_walks(path, next_x, next_y, length + 1)
                # Backtrack: remove the point to explore other possibilities.
                path.remove((next_x, next_y))

    # The walk must have at least one step.
    if n > 0:
        # Optimization: Exploit the lattice symmetry.
        # We calculate walks starting with a move to the right (1, 0)
        # and multiply the result by 4.
        # The path starts at (0,0) and takes one step to (1,0).
        initial_path = {(0, 0), (1, 0)}
        # Start the search from the end of the first step.
        find_walks(initial_path, 1, 0, 1)
        # The total number of walks is 4 times the count found for one direction.
        total_walks = count[0] * 4
    elif n == 0:
        # A 0-step walk is a single point, so there is 1 such walk.
        total_walks = 1
    else:
        total_walks = 0

    # Print the final result in the format "a(n) = result".
    print(f"a({n}) = {total_walks}")

if __name__ == "__main__":
    solve()