import math

def find_planar_distances():
    """
    Calculates the possible normalized distances (r <= 3) between hard spheres
    whose centers lie on a 2D plane, considering hexagonal and square packings.
    """
    
    # Use a set to store unique distances found.
    # Distances are normalized by the sphere diameter.
    distances = set()
    limit = 3.0
    
    # The search range for integer coefficients n and m.
    # A range of [-3, 3] is sufficient since for |n|>3 or |m|>3, r will be > 3.
    search_range = range(-limit, limit + 1)
    
    # --- Case 1: Hexagonal Lattice ---
    # The squared distance is r^2 = n^2 + nm + m^2.
    for n in search_range:
        for m in search_range:
            # Skip the origin (distance to itself)
            if n == 0 and m == 0:
                continue
            
            # Using r_sq for comparison is more robust against floating-point issues
            r_sq = float(n**2 + n * m + m**2)
            
            # Check if the distance is within the limit (with a small tolerance)
            if r_sq <= limit**2 + 1e-9:
                distances.add(math.sqrt(r_sq))

    # --- Case 2: Square Lattice ---
    # The squared distance is r^2 = n^2 + m^2.
    for n in search_range:
        for m in search_range:
            # Skip the origin
            if n == 0 and m == 0:
                continue
            
            r_sq = float(n**2 + m**2)
            
            # Check if the distance is within the limit (with a small tolerance)
            if r_sq <= limit**2 + 1e-9:
                distances.add(math.sqrt(r_sq))

    # Convert the set to a list and sort it
    sorted_distances = sorted(list(distances))
    
    # --- Output the final result ---
    print("The set of possible normalized distances r for spheres on a 2D plane (r <= 3) is:")
    
    # Format the output as an equation string with each number included
    equation_str = "r = " + ", ".join([f"{d:.2f}" for d in sorted_distances])
    print(equation_str)

# Execute the function
find_planar_distances()
<<<r = 1.00, 1.41, 1.73, 2.00, 2.24, 2.65, 2.83, 3.00>>>