import math

def find_planar_sphere_distances():
    """
    Calculates and prints the possible normalized distances (r <= 3) for jammed
    hard spheres whose centers are constrained to a 2D plane.

    The calculation is based on the geometry of a 2D hexagonal lattice, which represents
    the densest packing. The position of any sphere relative to a central sphere
    at the origin can be described by integer combinations (n, m) of two basis
    vectors. The squared distance r^2 is given by the formula r^2 = n^2 + m^2 + n*m,
    where the distance is normalized by the sphere diameter.
    """
    # The maximum normalized distance of interest.
    r_max = 3.0
    # Use a set to store the unique squared distances to avoid duplicates.
    r_squared_values = set()

    # The search range for integers n and m. A range up to r_max is sufficient.
    search_limit = int(r_max)

    for n in range(-search_limit, search_limit + 1):
        for m in range(-search_limit, search_limit + 1):
            # Skip the origin point (distance to itself).
            if n == 0 and m == 0:
                continue

            # Calculate the squared distance using the hexagonal lattice formula.
            r_squared = float(n*n + m*m + n*m)

            # Add the value if it's within the desired range (r <= 3 => r^2 <= 9).
            if r_squared <= r_max**2:
                r_squared_values.add(r_squared)

    # Convert the set of squared values to a sorted list of distances.
    distances = sorted([math.sqrt(r2) for r2 in r_squared_values])

    print("The set of possible distances (r) for r <= 3.00 are:")
    # Print each distance as a simple equation, formatted to two decimals.
    for r in distances:
        print(f"r = {r:.2f}")

find_planar_sphere_distances()
<<<[1.00, 1.73, 2.00, 2.65, 3.00]>>>