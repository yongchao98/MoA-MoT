import math

def find_hard_sphere_distances():
    """
    Calculates the possible distances between hard spheres on a 2D plane for r <= 3.

    This function models the sphere arrangement as a hexagonal lattice, which is the
    densest packing in 2D. The distances correspond to the distances between lattice points.
    A lattice point's position relative to a central particle at the origin can be
    represented by integer coefficients (n, m). The squared distance r^2 is given by
    the formula n^2 + nm + m^2.

    The function finds all unique distances r for which 0 < r <= 3, where r is
    normalized by the sphere diameter.
    """

    # Use a set to store unique distances to avoid duplicates.
    distances = set()
    limit = 3
    r_sq_limit = limit * limit

    # Iterate through a range of integer coefficients n and m. A range up to
    # the limit is sufficient to find all distances within the boundary.
    for n in range(-limit, limit + 1):
        for m in range(-limit, limit + 1):
            # Skip the case of the central particle itself (distance is 0).
            if n == 0 and m == 0:
                continue

            # Calculate the squared distance for a hexagonal lattice.
            r_sq = float(n**2 + n * m + m**2)

            # Add the distance to the set if it's within the limit r <= 3.
            if r_sq <= r_sq_limit:
                r = math.sqrt(r_sq)
                # To handle floating point inaccuracies, we round the result before adding
                # to the set. This ensures that distances like sqrt(3) from (1,1) and (2,-1)
                # are treated as the same.
                distances.add(round(r, 5))

    # Convert the set to a list and sort it for a clean, ordered output.
    sorted_distances = sorted(list(distances))

    # Print the final list of distances, formatted to two decimal places.
    print("The set of possible distances r (normalized by sphere diameter d) for r <= 3 are:")
    for r in sorted_distances:
        print(f"{r:.2f}")

if __name__ == '__main__':
    find_hard_sphere_distances()