def is_non_planar(n, m):
    """
    Checks if an n x m board generates a non-planar super-knight graph.
    Based on the known result, this is true if and only if n >= 5 and m >= 5.
    Args:
        n (int): number of rows, must be >= 4
        m (int): number of columns, must be >= 4
    Returns:
        bool: True if the graph is non-planar, False otherwise.
    """
    return n >= 5 and m >= 5

def find_supremum_size():
    """
    Finds the supremum of sizes nm for which all smaller or equal-sized boards are planar.
    This is equivalent to finding the minimum area of a non-planar board.
    """
    min_non_planar_area = -1
    
    # We search for the smallest area k >= 16 that can be formed by a non-planar rectangle
    # The smallest possible rectangle is 4x4=16
    k = 16 
    while True:
        # To check area k, we must check all its factor pairs (n,m)
        # where n, m >= 4.
        found_non_planar_for_k = False
        # Iterate through possible n values up to sqrt(k)
        for n in range(4, int(k**0.5) + 1):
            if k % n == 0:
                m = k // n
                # Check if this n x m rectangle is non-planar
                if is_non_planar(n, m):
                    min_non_planar_area = k
                    found_non_planar_for_k = True
                    break
        
        if found_non_planar_for_k:
            break
        k += 1
        
    return min_non_planar_area

# The first non-planar configuration occurs at n=5, m=5.
min_n = 5
min_m = 5
supremum = find_supremum_size()

print(f"A super-knight graph on an n x m board is non-planar if and only if n >= {min_n} and m >= {min_m}.")
print(f"The smallest area for a non-planar graph is therefore for a {min_n} x {min_m} board.")
print(f"Smallest non-planar area = {min_n} * {min_m} = {supremum}")
print(f"Any board with area strictly less than {supremum} is planar. The set of all k such that any board with area <= k is planar is (-inf, {supremum}).")
print(f"The supremum of this set is {supremum}.")
