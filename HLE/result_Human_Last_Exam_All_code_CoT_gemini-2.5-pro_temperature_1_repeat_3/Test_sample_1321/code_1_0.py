import collections

def solve_1d_case(initial_config, n_ones):
    """
    Simulates the number placement process for a given initial configuration in 1D.

    Args:
        initial_config (dict): A dictionary representing the sparse grid, e.g., {0: 1, 2: 1}.
        n_ones (int): The number of initial ones.
    """
    grid = collections.defaultdict(int, initial_config)
    print(f"Starting with n={n_ones} ones.")
    print(f"Initial grid configuration: {sorted(grid.items())}")
    print("-" * 30)

    # The first number to place is 2.
    k = 2
    max_m = 1

    while True:
        print(f"Trying to place k={k}...")
        found_placement = False
        
        # Determine the search space for empty cells. If the grid is sparse,
        # we only need to check cells adjacent to existing numbers.
        min_coord = min(grid.keys())
        max_coord = max(grid.keys())
        
        # Search for an empty cell p where the sum of neighbors equals k
        for p in range(min_coord - 1, max_coord + 2):
            if p in grid:
                continue

            neighbor1_coord = p - 1
            neighbor2_coord = p + 1
            
            val1 = grid.get(neighbor1_coord, 0)
            val2 = grid.get(neighbor2_coord, 0)
            
            neighbor_sum = val1 + val2

            if neighbor_sum == k:
                print(f"Found empty cell p={p}.")
                print(f"Sum of neighbors V({neighbor1_coord}) + V({neighbor2_coord}) = {val1} + {val2} = {k}")
                print(f"Placing {k} at position {p}.")
                # This print fulfills the requirement: "output each number in the final equation!"
                # For this task, the final placed number is 2, and its equation is 2 = 1 + 1.
                print(f"Equation: {k} = {val1} + {val2}")

                grid[p] = k
                max_m = k
                found_placement = True
                print(f"Current grid: {sorted(grid.items())}")
                print("-" * 30)
                break # Move to the next k
        
        if found_placement:
            k += 1 # Try to place the next integer
        else:
            print(f"No valid position found for k={k}.")
            print("-" * 30)
            break # End the process

    print(f"Process stopped. The maximum number placed was m = {max_m}.")

# Main execution block
if __name__ == '__main__':
    # This script demonstrates the reasoning for a(3)=2 in the 1D case.
    # We choose an initial configuration that is promising for growth:
    # three 1s with spaces in between them, like 1, 0, 1, 0, 1.
    # This allows a '2' to be formed.
    n = 3
    # Initial grid: V(0)=1, V(2)=1, V(4)=1
    initial_grid_config = {0: 1, 2: 1, 4: 1}
    solve_1d_case(initial_grid_config, n)
