import sys

def solve():
    """
    Calculates the minimum number of initially infected sites to infect a 14x14 grid.
    The method is based on the necessary conditions for bootstrap percolation with threshold 3.
    """
    n = 14

    # Condition 1: The initial set must hit every 2x2 square.
    # We use a pattern known to satisfy this: infect sites (r, c) where (r - c) is a multiple of 3.
    base_pattern_sites = set()
    for r in range(n):
        for c in range(n):
            if (r - c) % 3 == 0:
                base_pattern_sites.add((r, c))

    base_count = len(base_pattern_sites)

    # Condition 2: The four corners must be infected.
    corners = {
        (0, 0),
        (0, n - 1),
        (n - 1, 0),
        (n - 1, n - 1)
    }

    # Combine the sets. The union automatically handles corners already in the base pattern.
    final_set = base_pattern_sites.union(corners)
    final_count = len(final_set)
    
    # The number of corners we needed to add is the difference between the final count
    # and the base pattern count.
    added_corners_count = final_count - base_count
    
    # Output the result following the required format.
    print(f"The calculation is based on ensuring all 2x2 squares and all four corners are infected.")
    print(f"Number of sites from the pattern (r-c) % 3 == 0: {base_count}")
    print(f"Number of corners that must be added to the pattern: {added_corners_count}")
    print(f"The final equation for the total number of sites is:")
    print(f"{base_count} + {added_corners_count} = {final_count}")

solve()
<<<68>>>