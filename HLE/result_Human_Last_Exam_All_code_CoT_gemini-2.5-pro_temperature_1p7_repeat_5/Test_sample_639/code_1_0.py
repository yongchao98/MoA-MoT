import math

def find_packing_distances():
    """
    Calculates the possible normalized distances between hard spheres on a 2D plane
    for r <= 3, based on hexagonal and square lattice structures. These correspond
    to the expected peaks in the radial distribution function g(r).
    """
    
    # Use a set to automatically store unique distance values.
    distances_set = set()
    limit = 3.0
    
    # Define a search range for the integer lattice vectors (n, m).
    # A range from -3 to 3 is sufficient, as combinations with larger integers
    # will result in distances r > 3.
    int_range = range(-3, 4)
    
    for n in int_range:
        for m in int_range:
            # Skip the origin point (distance to self).
            if n == 0 and m == 0:
                continue
                
            # --- Case 1: Hexagonal lattice ---
            # The squared normalized distance is r^2 = n^2 + m^2 + nm.
            r_sq_hex = float(n**2 + m**2 + n * m)
            if r_sq_hex <= limit**2:
                distances_set.add(math.sqrt(r_sq_hex))
                
            # --- Case 2: Square lattice ---
            # The squared normalized distance is r^2 = n^2 + m^2.
            r_sq_square = float(n**2 + m**2)
            if r_sq_square <= limit**2:
                distances_set.add(math.sqrt(r_sq_square))

    # Convert the set to a list and sort it for ordered output.
    sorted_distances = sorted(list(distances_set))
    
    print("The set of possible normalized distances r <= 3 are:")
    
    for r in sorted_distances:
        # We check if r^2 is close to an integer to find the original equation.
        # A small tolerance is used to handle floating point inaccuracies.
        r_sq = r**2
        r_sq_int = round(r_sq)
        
        # If r is a whole number (e.g., sqrt(1), sqrt(4), sqrt(9)).
        if abs(r - round(r)) < 1e-9:
            print(f"r = {r:.2f}")
        # If r is the square root of an integer.
        elif abs(r_sq - r_sq_int) < 1e-9:
            # Provide a more descriptive form for common radicals like sqrt(8).
            if r_sq_int == 8:
                print(f"r = sqrt(8) = 2*sqrt(2) = {r:.2f}")
            else:
                print(f"r = sqrt({r_sq_int}) = {r:.2f}")
        # This case should not be reached with the current models.
        else:
            print(f"r = {r:.2f}")

# Execute the function to print the results.
find_packing_distances()