import math

def find_packing_distances():
    """
    Calculates the set of possible normalized distances (r) between hard spheres
    on a 2D plane for r <= 3.

    This is done by considering the distances found in the two most common 2D
    lattices: square and hexagonal. A jammed amorphous system will exhibit local
    structures resembling both. The distances are normalized by the sphere diameter.
    """
    
    # Use a set to automatically handle duplicate distances
    distances = set()

    # Define a search range for lattice coordinates (n, m).
    # A range up to 3 is sufficient since r <= 3.
    search_range = range(-3, 4)

    # --- Hexagonal Lattice Calculations ---
    # The squared distance from the origin to a point (n, m) in a hexagonal
    # lattice (with side length 1) is given by r^2 = n^2 + nm + m^2.
    for n in search_range:
        for m in search_range:
            # Skip the origin particle
            if n == 0 and m == 0:
                continue
            
            r_squared = float(n**2 + n * m + m**2)
            r = math.sqrt(r_squared)
            
            # Add distance to the set if it's within the desired range (r<=3)
            # Use a small tolerance for floating point comparison
            if r <= 3.00001:
                distances.add(r)

    # --- Square Lattice Calculations ---
    # The squared distance from the origin to a point (n, m) in a square
    # lattice (with side length 1) is given by r^2 = n^2 + m^2.
    for n in search_range:
        for m in search_range:
            # Skip the origin particle
            if n == 0 and m == 0:
                continue
            
            r_squared = float(n**2 + m**2)
            r = math.sqrt(r_squared)
            
            # Add distance to the set if it's within the desired range (r<=3)
            if r <= 3.00001:
                distances.add(r)
    
    # Convert the set to a sorted list for ordered output
    sorted_distances = sorted(list(distances))
    
    print("The set of possible normalized distances r <= 3 are:")
    # Print each unique distance formatted to two decimal places
    for r in sorted_distances:
        print(f"{r:.2f}")

if __name__ == "__main__":
    find_packing_distances()