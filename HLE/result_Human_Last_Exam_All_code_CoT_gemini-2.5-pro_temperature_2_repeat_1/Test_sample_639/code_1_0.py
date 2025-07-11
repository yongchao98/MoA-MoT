import math

def find_hs_distances():
    """
    Calculates the set of possible normalized distances (r <= 3) between
    two hard spheres in a 2D jammed packing.

    The calculation is based on the distances in a hexagonal lattice, which is the
    densest packing in 2D. The squared distance r^2 from a central sphere to any
    other sphere can be described by r^2 = m^2 + m*n + n^2, where m and n are integers.
    """
    # The maximum normalized distance to consider
    limit = 3

    # Use a set to store unique squared distances to avoid duplicates
    squared_distances = set()

    # Iterate through a grid of integer indices (m, n) to find all possible
    # sphere positions. A range up to 'limit' is sufficient to find all
    # distances r <= limit.
    for m in range(-limit, limit + 1):
        for n in range(-limit, limit + 1):
            # Skip the case of the central particle itself (distance zero)
            if m == 0 and n == 0:
                continue

            # Calculate the squared distance based on the hexagonal lattice formula
            r_sq = m**2 + m * n + n**2

            # Add the squared distance to our set if it's within the limit
            if r_sq <= limit**2:
                squared_distances.add(r_sq)

    # Convert the set of squared distances to a sorted list of actual distances
    distances = sorted([math.sqrt(r_sq) for r_sq in squared_distances])

    # Print the final result in the required format
    print("The set of possible normalized distances r <= 3, representing the main peaks in g(r), is:")
    
    # We create a formatted string for the final output, as if it were part of an equation or result list.
    # The instructions state: "Remember in the final code you still need to output each number in the final equation!"
    result_list = []
    for r in distances:
        result_list.append(f"{r:.2f}")

    print("r = {" + ", ".join(result_list) + "}")

if __name__ == '__main__':
    find_hs_distances()