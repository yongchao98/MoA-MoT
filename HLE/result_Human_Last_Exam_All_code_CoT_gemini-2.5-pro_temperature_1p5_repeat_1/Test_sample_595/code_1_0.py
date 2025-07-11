import math

def calculate_k(i, j):
    """Calculates the number of squares k based on integer parts i and j."""
    return 3 * (i + j) + abs(j - i) + 1

def solve():
    """
    Finds the maximum number of squares the triangle's perimeter can pass through.
    """
    max_k = 0
    best_pair = (0, 0)

    # i = floor(18 * cos(phi))
    # j = floor(18 * sin(phi))
    # This implies i^2 + j^2 must be less than 18^2 = 324.
    # We iterate through all possible integer values for i.
    for i in range(18):
        # For a given i, find the maximum possible j.
        # i^2 + j^2 < 324 => j^2 < 324 - i^2
        # j < sqrt(324 - i^2)
        j_limit = math.floor(math.sqrt(324 - 1 - i*i)) # -1 for safety margin as it's a strict inequality
        for j in range(j_limit + 1):
            k = calculate_k(i, j)
            if k > max_k:
                max_k = k
                best_pair = (i, j)

    i, j = best_pair
    # The axis-aligned case corresponds to i=18, j=0 or i=0, j=18, which our loop misses.
    # Let's check that case manually. Vertices (e,e), (18+e,e), (e,18+e)
    # Side 1 (leg): dx=18, dy=0. Crossings = 18.
    # Side 2 (leg): dx=0, dy=18. Crossings = 18.
    # Side 3 (hypotenuse): dx=-18, dy=18. Crossings = 18+18 = 36.
    # Total k_aligned = 18 + 18 + 36 = 72.
    k_aligned = 72
    print(f"The number of squares crossed in the axis-aligned case is 18 + 18 + 36 = {k_aligned}.")
    if k_aligned > max_k:
        max_k = k_aligned
        best_pair = (18, 0)
    
    print(f"The optimal integer pair (i, j) found is {best_pair} (or {(best_pair[1], best_pair[0])}).")
    print(f"The calculation for k with this pair is: 3*({i}+{j}) + abs({j}-{i}) + 1 = {max_k}.")
    print(f"The largest number of coordinate grid squares is {max_k}.")


solve()