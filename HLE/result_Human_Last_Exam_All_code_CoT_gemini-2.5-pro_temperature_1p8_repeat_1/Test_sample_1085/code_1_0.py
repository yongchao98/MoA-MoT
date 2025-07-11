import collections

# Using a larger set of primes to ensure correctness for any intermediate steps.
PRIMES_AND_ONE = [1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

# Create mappings for efficient lookup of prime indices
PRIME_TO_INDEX = {p: i for i, p in enumerate(PRIMES_AND_ONE)}
INDEX_TO_PRIME = {i: p for i, p in enumerate(PRIMES_AND_ONE)}

# Memoization cache for the recursive path counting
memo = {}

def get_adjacent_nodes(node):
    """
    Finds all adjacent nodes based on the "exactly one prime in between" rule.
    This is equivalent to the prime index differing by 2.
    """
    px, py = node
    if px not in PRIME_TO_INDEX or py not in PRIME_TO_INDEX:
        return []

    idx_x = PRIME_TO_INDEX[px]
    idx_y = PRIME_TO_INDEX[py]
    
    adjacent = []

    # Horizontal neighbors
    for diff in [-2, 2]:
        next_idx_x = idx_x + diff
        if next_idx_x in INDEX_TO_PRIME:
            adjacent.append((INDEX_TO_PRIME[next_idx_x], py))

    # Vertical neighbors
    for diff in [-2, 2]:
        next_idx_y = idx_y + diff
        if next_idx_y in INDEX_TO_PRIME:
            adjacent.append((px, INDEX_TO_PRIME[next_idx_y]))
            
    return adjacent

def find_paths_recursive(start_node, moves_left, target_node, current_path):
    """
    Recursively finds and returns all valid paths.
    """
    if moves_left == 0:
        if start_node == target_node:
            return [current_path]
        else:
            return []

    all_paths = []
    neighbors = get_adjacent_nodes(start_node)
    
    for neighbor in neighbors:
        # To find distinct paths, we can revisit nodes.
        found_paths = find_paths_recursive(neighbor, moves_left - 1, target_node, current_path + [neighbor])
        all_paths.extend(found_paths)
        
    return all_paths

def solve():
    """
    Main function to find and print the number of distinct prime paths.
    """
    start_node = (1, 1)
    target_node = (5, 7)
    total_moves = 4
    
    print(f"Finding all distinct paths from {start_node} to {target_node} in {total_moves} moves.\n")

    # We use a list to keep track of the path
    initial_path = [start_node]
    
    all_paths = find_paths_recursive(start_node, total_moves, target_node, initial_path)
    
    path_count = len(all_paths)
    
    print(f"Found {path_count} distinct paths.")
    
    if path_count > 0:
        print("\nThe paths are:")
        for path in all_paths:
            # We want to print each number in the equation.
            path_str = " -> ".join(map(str, path))
            print(path_str)
    
    total = 0
    final_paths = []
    
    if path_count > 0:
        for path in all_paths:
            final_paths.append("1")
        total = "+".join(final_paths)
        print(f"\nThe final equation for the number of paths is: {total} = {path_count}")
    else:
        print(f"\nThe final equation for the number of paths is: 0 = 0")


if __name__ == "__main__":
    solve()