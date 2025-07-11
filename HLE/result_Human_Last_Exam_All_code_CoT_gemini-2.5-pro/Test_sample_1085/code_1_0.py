import collections

# Using a dictionary for memoization to store results of (node, moves_left)
memo = {}
# The ordered list of allowed coordinates defines the grid's structure
ALLOWED_COORDS = [1, 2, 3, 5, 7]
# A mapping from coordinate value to its index for quick lookups
COORD_TO_INDEX = {val: i for i, val in enumerate(ALLOWED_COORDS)}

def get_neighbors(node):
    """
    Calculates the valid neighbors of a node in the PrimeGrid+1.
    A neighbor is one step away horizontally or vertically on the prime-indexed grid.
    """
    x, y = node
    neighbors = []
    x_idx = COORD_TO_INDEX.get(x)
    y_idx = COORD_TO_INDEX.get(y)

    # Horizontal neighbors
    if x_idx > 0:
        neighbors.append((ALLOWED_COORDS[x_idx - 1], y))
    if x_idx < len(ALLOWED_COORDS) - 1:
        neighbors.append((ALLOWED_COORDS[x_idx + 1], y))
    
    # Vertical neighbors
    if y_idx > 0:
        neighbors.append((x, ALLOWED_COORDS[y_idx - 1]))
    if y_idx < len(ALLOWED_COORDS) - 1:
        neighbors.append((x, ALLOWED_COORDS[y_idx + 1]))
        
    return neighbors

def count_paths(current, target, moves_left):
    """
    Recursively counts the number of paths from the current node to the target
    with a specific number of moves left.
    """
    # Check memoization table first
    if (current, moves_left) in memo:
        return memo[(current, moves_left)]

    # Base case: if no moves are left
    if moves_left == 0:
        # Return 1 if we are at the target, 0 otherwise
        return 1 if current == target else 0

    # Recursive step: sum the paths from all neighbors
    path_count = 0
    for neighbor in get_neighbors(current):
        path_count += count_paths(neighbor, target, moves_left - 1)
    
    # Store the result in the memoization table before returning
    memo[(current, moves_left)] = path_count
    return path_count

def solve_prime_path_problem():
    """
    Solves the main problem by setting up the initial call and printing the results.
    """
    start_node = (1, 1)
    end_node = (5, 7)
    total_moves = 4
    
    print(f"Finding paths from {start_node} to {end_node} in {total_moves} moves.")
    print("-" * 30)

    # A path must start by moving to a neighbor of the start_node.
    # The total number of paths is the sum of paths starting from each neighbor.
    initial_neighbors = get_neighbors(start_node)
    
    # The first move is from (1,1) to either (1,2) or (2,1).
    # We need to find the number of paths of length 3 from these neighbors.
    moves_remaining = total_moves - 1
    
    count_from_1_2 = count_paths((1, 2), end_node, moves_remaining)
    count_from_2_1 = count_paths((2, 1), end_node, moves_remaining)
    
    total_paths = count_from_1_2 + count_from_2_1

    print(f"Number of paths starting with (1,1) -> (1,2): {count_from_1_2}")
    print(f"Number of paths starting with (1,1) -> (2,1): {count_from_2_1}")
    print("\nFinal equation for the total number of paths:")
    print(f"Total Paths = {count_from_1_2} + {count_from_2_1} = {total_paths}")

solve_prime_path_problem()
<<<0>>>