import math

def find_planar_distances():
    """
    Calculates the possible distances between centers of hard spheres
    arranged on a single 2D plane in a jammed configuration, for r <= 3.

    This is equivalent to finding the lengths of vectors on a hexagonal lattice.
    The squared distance r^2 (normalized by diameter d) between any two lattice
    points is given by r^2 = n^2 + nm + m^2, where n and m are integers.
    """
    # The maximum distance to consider is r=3, so r^2=9
    r_max_sq = 3.0 * 3.0

    # Use a set to store unique distances
    distances = set()

    # Define a search range for the integer indices n and m.
    # A small integer range is sufficient to find all distances up to r=3.
    search_range = 3

    # Iterate through lattice indices n and m
    for n in range(-search_range, search_range + 1):
        for m in range(-search_range, search_range + 1):
            # Skip the (0,0) case which corresponds to a zero distance
            if n == 0 and m == 0:
                continue

            # Calculate the squared distance
            r_sq = float(n*n + n*m + m*m)

            # Check if the distance is within the specified limit
            if r_sq <= r_max_sq:
                dist = math.sqrt(r_sq)
                # We round here to handle potential floating point inaccuracies
                # before adding to the set, ensuring uniqueness.
                distances.add(round(dist, 5))

    # Sort the collected distances
    sorted_distances = sorted(list(distances))

    # Print each distance formatted to two decimal places
    print("The set of normalized distances r <= 3 are:")
    for r in sorted_distances:
        print(f"{r:.2f}")

if __name__ == "__main__":
    find_planar_distances()