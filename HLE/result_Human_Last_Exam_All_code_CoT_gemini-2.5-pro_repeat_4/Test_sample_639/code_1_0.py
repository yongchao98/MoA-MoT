import math

def find_planar_distances():
    """
    Calculates the possible normalized distances (r) between hard spheres
    on a 2D hexagonal lattice for r <= 3.

    In a 2D hexagonal lattice, the squared distance from a central point
    to any other point (n, m) is given by r^2 = n^2 + nm + m^2.
    """
    # Use a set to store unique squared distances to avoid duplicates
    r_squared_values = set()
    
    # The range for n and m can be limited, as r^2 must be <= 9.
    # A range of -3 to 3 for both n and m is sufficient.
    limit = 3
    for n in range(-limit, limit + 1):
        for m in range(-limit, limit + 1):
            # Exclude the origin point (0,0)
            if n == 0 and m == 0:
                continue
            
            r_sq = n**2 + n * m + m**2
            
            # We are interested in distances r <= 3, so r^2 <= 9
            if r_sq <= 9:
                r_squared_values.add(r_sq)
                
    # Sort the unique squared distances
    sorted_r_sq = sorted(list(r_squared_values))
    
    print("The set of possible distances r for r <= 3, based on a 2D hexagonal packing, are:")
    
    final_distances = []
    # Calculate and print the final r values
    for r_sq in sorted_r_sq:
        r = math.sqrt(r_sq)
        final_distances.append(f"{r:.2f}")
        # The output format shows each number in the final equation as requested
        print(f"r = sqrt({r_sq}) = {r:.2f}")
    
    # Return the final values for the answer block
    return ", ".join(final_distances)

if __name__ == '__main__':
    find_planar_distances()
