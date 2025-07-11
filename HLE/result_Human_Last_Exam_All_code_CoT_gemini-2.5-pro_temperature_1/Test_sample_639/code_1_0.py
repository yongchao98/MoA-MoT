import math

def find_hs_distances():
    """
    Calculates the set of possible normalized distances (r <= 3) between
    hard spheres in a 2D jammed system.

    The calculation considers local sphere arrangements based on the square
    and hexagonal lattices, which represent the most common packing motifs.
    """
    # Use a set to store unique squared distances (r^2) to avoid duplicates.
    # We store r^2 to maintain precision before the final square root.
    r_squared_values = set()

    # Define the range for integers m and n. A range up to 3 is sufficient
    # because for r=3, the largest m or n could be is 3 (e.g., 3^2=9).
    integer_range = range(-3, 4)

    # --- Part 1: Hexagonal Lattice Calculations ---
    # The squared distance formula is r^2 = m^2 + n^2 + mn.
    for m in integer_range:
        for n in integer_range:
            # Skip the origin point (distance to self)
            if m == 0 and n == 0:
                continue
            
            r_sq = float(m**2 + n**2 + m * n)
            
            # We are interested in distances r <= 3, which means r^2 <= 9.
            if r_sq <= 9.000001: # Use a small tolerance for float comparison
                r_squared_values.add(r_sq)

    # --- Part 2: Square Lattice Calculations ---
    # The squared distance formula is r^2 = m^2 + n^2.
    for m in integer_range:
        for n in integer_range:
            # Skip the origin point
            if m == 0 and n == 0:
                continue
            
            r_sq = float(m**2 + n**2)
            
            # We are interested in distances r <= 3, so r^2 <= 9.
            if r_sq <= 9.000001:
                r_squared_values.add(r_sq)

    # Sort the unique squared distances in ascending order
    sorted_r_squared = sorted(list(r_squared_values))

    print("The set of possible normalized distances r for r <= 3 are:")
    
    # --- Part 3: Print the Final Results ---
    # We loop through the sorted list of r^2 values, calculate r,
    # and print the final formatted equation for each distance.
    for r_sq in sorted_r_squared:
        r = math.sqrt(r_sq)
        # Output each number in a final equation as requested.
        print(f"r = {r:.2f}")

# Execute the function
find_hs_distances()