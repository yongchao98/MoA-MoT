import math

def find_jammed_distances():
    """
    Calculates the possible normalized distances (r <= 3) between hard spheres
    packed on a 2D plane by considering hexagonal and square lattice arrangements.
    """

    # The maximum normalized distance is 3, so the maximum r^2 is 9.
    max_r_squared = 9.0

    # Using a set to store unique squared distances to avoid duplicates
    r_squared_values = set()

    # We iterate through possible integer lattice vector components (n, m).
    # A search range up to 3 for n and m is sufficient, since if n or m > 3,
    # n^2 or m^2 would be greater than 9 for most cases.
    # range(4) -> 0, 1, 2, 3
    limit = 4 

    for n in range(limit):
        for m in range(limit):
            # We skip the origin point (distance to self)
            if n == 0 and m == 0:
                continue

            # --- Case 1: Hexagonal Lattice ---
            # The squared distance for a hexagonal lattice can be calculated as r^2 = n^2 + nm + m^2
            # for integer vectors (n, m) in a hexagonal coordinate system.
            r_sq_hex = float(n**2 + n * m + m**2)
            if r_sq_hex <= max_r_squared:
                r_squared_values.add(r_sq_hex)

            # --- Case 2: Square Lattice ---
            # The squared distance for a square lattice is r^2 = n^2 + m^2.
            r_sq_square = float(n**2 + m**2)
            if r_sq_square <= max_r_squared:
                r_squared_values.add(r_sq_square)

    # Calculate the actual distances by taking the square root of the unique squared values
    distances = sorted([math.sqrt(r_sq) for r_sq in r_squared_values])

    # Print the final results
    print("The set of possible normalized distances r for particles on a plane with r <= 3.00 are:")
    for dist in distances:
        # The final requirement is to output each number from the final calculation.
        print(f"{dist:.2f}")

if __name__ == "__main__":
    find_jammed_distances()
