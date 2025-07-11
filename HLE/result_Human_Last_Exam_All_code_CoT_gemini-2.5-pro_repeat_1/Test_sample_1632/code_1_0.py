def count_walks(x, y, steps_left, visited):
    """
    Recursively counts the number of self-avoiding walks of a given length
    from a starting point (x, y), avoiding points in the 'visited' set.
    """
    # Base case: If there are no steps left, we have found one valid walk.
    if steps_left == 0:
        return 1

    count = 0
    # Define the four possible moves on a Manhattan lattice.
    moves = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

    # Explore each possible move.
    for next_x, next_y in moves:
        # If the next point has not been visited yet...
        if (next_x, next_y) not in visited:
            # ...add it to the visited set and recurse.
            visited.add((next_x, next_y))
            count += count_walks(next_x, next_y, steps_left - 1, visited)
            # Backtrack: remove the point to allow other paths to visit it.
            visited.remove((next_x, next_y))
            
    return count

def main():
    """
    Calculates a(10), the number of 10-step self-avoiding walks.
    """
    n = 10

    # Due to the symmetry of the lattice, we can calculate the number of walks
    # starting with a single direction (e.g., right) and multiply the result by 4.
    # The walk starts at (0,0). The first step is to (1,0).
    start_x, start_y = 1, 0
    steps_to_take = n - 1
    
    # The initial visited set contains the origin (0,0) and the first step's destination (1,0).
    initial_visited = {(0, 0), (start_x, start_y)}

    # Count the number of valid (n-1)-step walks from the first step's position.
    walks_for_one_direction = count_walks(start_x, start_y, steps_to_take, initial_visited)
    
    # The total number of walks is 4 times this count.
    total_walks = 4 * walks_for_one_direction
    
    # Print the final result in the format "a(n) = result".
    print(f"a({n}) = {total_walks}")

if __name__ == "__main__":
    main()