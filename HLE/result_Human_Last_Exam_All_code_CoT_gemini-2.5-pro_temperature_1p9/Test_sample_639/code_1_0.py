import math

def find_hs_distances():
    """
    Calculates the possible normalized distances between hard spheres
    in a 2D jammed packing for distances r <= 3.

    The sphere centers form a triangular lattice. The squared distance r^2
    from a central sphere to any other sphere at lattice coordinate (n, m)
    is given by r^2 = n^2 + nm + m^2.
    """
    # Use a set to automatically store only unique distances
    unique_distances = set()

    # Iterate through a range of integer lattice coordinates (n, m)
    # A range of -5 to 5 is more than sufficient to find all r <= 3.
    limit = 5
    for n in range(-limit, limit + 1):
        for m in range(-limit, limit + 1):
            # Skip the case of the distance to the sphere itself
            if n == 0 and m == 0:
                continue

            # Calculate the squared distance using the formula for a triangular lattice
            r_squared = n**2 + n * m + m**2

            # We are interested in distances r <= 3, which is r^2 <= 9
            # Add a small tolerance for floating point comparisons
            if r_squared <= 9.0000001:
                distance = math.sqrt(r_squared)
                unique_distances.add(distance)

    # Sort the final list of unique distances
    sorted_distances = sorted(list(unique_distances))

    print("The possible normalized distances r for r <= 3 are:")
    for dist in sorted_distances:
        print(f"{dist:.2f}")

if __name__ == '__main__':
    find_hs_distances()