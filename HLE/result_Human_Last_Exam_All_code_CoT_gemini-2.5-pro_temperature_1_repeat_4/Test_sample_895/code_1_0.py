import math

def calculate_max_hausdorff_distance(a, n):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B.

    Args:
        a (list of float): The edge lengths of the polygon B.
        n (int): The number of vertices/edges of the polygon B.
    """
    if len(a) != n:
        print("Error: The number of edge lengths must be equal to n.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    max_params = {}

    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0:
            # This can happen with floating point inaccuracies for b_i near 0
            b_i_sq = 0
        b_i = math.sqrt(b_i_sq)

        if b_i == 0:
            # This case corresponds to a degenerate polygon, skip
            continue
            
        dist = (a_i * a_i_plus_1 * sin_phi) / (2 * b_i)

        if dist > max_dist:
            max_dist = dist
            max_params = {
                'i': i,
                'a_i': a_i,
                'a_i_plus_1': a_i_plus_1,
                'phi': phi,
                'cos_phi': cos_phi,
                'sin_phi': sin_phi,
                'b_i': b_i
            }
            
    if max_dist == -1:
        print("Could not compute the distance. The polygon might be degenerate.")
        return
        
    print("The formula for the largest possible Hausdorff distance is the maximum of:")
    print("H(i) = (a_i * a_{i+1} * sin(phi)) / (2 * b_i)")
    print("where b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)) and phi = 2*pi/n\n")

    i = max_params['i']
    a_i = max_params['a_i']
    a_i_plus_1 = max_params['a_i_plus_1']
    phi_val = max_params['phi']
    b_i_val = max_params['b_i']

    print(f"The maximum value is found for the vertex between edge {i} and edge {i+1}.")
    print("\n--- Numbers for the final equation ---")
    print(f"i = {i}")
    print(f"a_i = {a_i}")
    print(f"a_{i+1} = {a_i_plus_1}")
    print(f"n = {n}")
    print(f"phi (rad) = {phi_val:.4f}")
    print(f"sin(phi) = {max_params['sin_phi']:.4f}")
    print(f"cos(phi) = {max_params['cos_phi']:.4f}")
    print(f"b_i = sqrt({a_i}^2 + {a_i_plus_1}^2 + 2*{a_i}*{a_i_plus_1}*cos({phi_val:.4f})) = {b_i_val:.4f}")
    
    print("\n--- Final Calculation ---")
    print(f"Distance = ({a_i} * {a_i_plus_1} * {max_params['sin_phi']:.4f}) / (2 * {b_i_val:.4f})")
    print(f"Largest possible Hausdorff distance = {max_dist:.4f}")

if __name__ == '__main__':
    # Example usage: A rectangle with sides 3 and 4
    # For a rectangle, n=4. The normal vectors are at 0, 90, 180, 270 degrees.
    n_vertices = 4
    edge_lengths = [3.0, 4.0, 3.0, 4.0]
    
    # You can change the values here to test other polygons
    # Example: A regular hexagon with side length 5
    # n_vertices = 6
    # edge_lengths = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]

    calculate_max_hausdorff_distance(edge_lengths, n_vertices)
