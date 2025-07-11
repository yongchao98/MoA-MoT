import math

def find_packed_sphere_distances():
    """
    Calculates the set of possible normalized distances (r) between hard spheres
    in a 2D plane for r <= 3.

    This is done by considering the two primary regular packing lattices in 2D:
    square and hexagonal. A jammed system can have local arrangements
    corresponding to either geometry.
    """
    
    # Use a set to store unique distances to avoid duplicates
    distances = set()
    
    # The maximum integer coordinate to check. If r <= 3, r^2 <= 9.
    # For both n^2+m^2 and n^2+m^2+nm, coordinates will not exceed 3.
    # A search range of -3 to 3 is sufficient.
    coord_range = range(-3, 4)
    
    # Iterate through integer combinations (n, m) for the lattice vectors
    for n in coord_range:
        for m in coord_range:
            # Skip the origin particle (distance to itself is 0)
            if n == 0 and m == 0:
                continue
            
            # --- Case 1: Square lattice ---
            # The squared distance is r^2 = n^2 + m^2 (in units of diameter d)
            r_squared_square = float(n**2 + m**2)
            if r_squared_square <= 9:
                distances.add(math.sqrt(r_squared_square))

            # --- Case 2: Hexagonal lattice ---
            # The squared distance is r^2 = n^2 + m^2 + nm
            r_squared_hex = float(n**2 + m**2 + n * m)
            if r_squared_hex <= 9:
                distances.add(math.sqrt(r_squared_hex))

    # Sort the collected distances in ascending order
    sorted_distances = sorted(list(distances))
    
    print("The set of possible normalized distances r <= 3 are:")
    # Print each distance formatted to two decimal places
    for r in sorted_distances:
        # The final output needs to be the numbers of the final equation.
        # Here we print each component of the set.
        print(f"{r:.2f}")

if __name__ == "__main__":
    find_packed_sphere_distances()
