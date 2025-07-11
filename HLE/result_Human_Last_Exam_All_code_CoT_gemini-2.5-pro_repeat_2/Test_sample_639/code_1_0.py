import math

def find_sphere_distances():
    """
    Calculates the possible normalized distances between hard spheres on a 2D plane
    for distances r <= 3.

    The calculation is based on the geometry of a hexagonal lattice, which represents
    the densest packing of spheres in a plane. The squared normalized distance r^2
    from a central sphere to any other sphere is given by r^2 = n^2 + n*m + m^2,
    where n and m are integer coordinates.
    """

    # We are looking for distances r <= 3, so the maximum squared distance is 3*3 = 9.
    r_max_sq = 9

    # Use a set to automatically store only the unique squared distances.
    unique_r_sq_values = set()

    # Iterate through integer coordinates n and m. A range from -3 to 3 is
    # sufficient, as larger values would result in r^2 > 9.
    for n in range(-3, 4):
        for m in range(-3, 4):
            # Skip the case (0,0) which is the reference sphere itself.
            if n == 0 and m == 0:
                continue

            # Calculate the squared distance using the hexagonal lattice formula.
            r_sq = n**2 + n * m + m**2

            # If the distance is within our limit, add it to the set.
            if r_sq <= r_max_sq:
                unique_r_sq_values.add(r_sq)

    # Sort the unique squared distances for an ordered output.
    sorted_r_sq = sorted(list(unique_r_sq_values))

    print("The possible normalized distances r for r <= 3 are:")
    
    final_distances = []
    # Calculate and print each distance r from its squared value r_sq.
    for r_sq in sorted_r_sq:
        r = math.sqrt(r_sq)
        print(f"For r^2 = {r_sq}, the distance is r = sqrt({r_sq}) = {r:.2f}")
        final_distances.append(f"{r:.2f}")
        
    # This part is for the final answer format, not part of the primary output.
    # print(f"\nFinal set of distances: {{{', '.join(final_distances)}}}")

if __name__ == "__main__":
    find_sphere_distances()
