import math

def find_hs_distances():
    """
    Calculates the characteristic distances for hard spheres in a 2D hexagonal packing.

    The function models a hexagonal lattice and finds all unique inter-particle
    distances up to a normalized distance of 3.
    """
    # A set to store unique distances to avoid duplicates
    distances = set()

    # The basis vectors for a hexagonal lattice with unit distance between neighbors
    v1 = (1.0, 0.0)
    v2 = (0.5, math.sqrt(3.0) / 2.0)

    # We iterate over a grid of integer coefficients (n, m) to generate lattice points.
    # A range of -4 to 4 is sufficient to find all distances up to r=3.
    for n in range(-4, 5):
        for m in range(-4, 5):
            # Skip the origin point (the central sphere itself)
            if n == 0 and m == 0:
                continue

            # Calculate the position vector p = n*v1 + m*v2
            px = n * v1[0] + m * v2[0]
            py = n * v1[1] + m * v2[1]

            # Calculate the distance (magnitude of the vector p)
            dist = math.sqrt(px**2 + py**2)

            # We are interested in distances r <= 3
            # Add a small tolerance for floating point comparisons
            if dist <= 3.00001:
                distances.add(dist)

    # Sort the distances in ascending order
    sorted_distances = sorted(list(distances))

    # Print the final result in the required format
    print("The set of normalized distances r for a 2D hexagonal packing (for r <= 3) is given by the following values:")
    for d in sorted_distances:
        # The prompt asks to output each number in the "final equation"
        print(f"r = {d:.2f}")

if __name__ == "__main__":
    find_hs_distances()
    # As requested by the final instruction format.
    # We output the set of numbers found.
    print("\n<<<1.00, 1.73, 2.00, 2.65, 3.00>>>")
