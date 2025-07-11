import math

def calculate_max_hausdorff_distance(edge_lengths):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B.

    Args:
        edge_lengths (list of float): A list of the edge lengths a_1, ..., a_n of the polygon B.

    Returns:
        float: The largest possible Hausdorff distance.
    """
    n = len(edge_lengths)
    if n < 3:
        raise ValueError("A polygon must have at least 3 edges.")

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = 0.0

    print("Step-by-step calculation:")
    print(f"Number of sides n = {n}")
    print(f"Angle phi = 2*pi/n = {phi:.4f} radians")
    print("-" * 30)

    for i in range(n):
        a_i = edge_lengths[i]
        # The index (i + 1) % n ensures a wrap-around for the last edge
        a_i_plus_1 = edge_lengths[(i + 1) % n]

        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        b_i = math.sqrt(b_i_squared)

        # The maximal distance for the i-th vertex is h_i = (a_i * a_{i+1} * sin(phi)) / (2 * b_i)
        numerator = a_i * a_i_plus_1 * sin_phi
        h_i = numerator / (2 * b_i)

        print(f"For i = {i+1}:")
        print(f"  a_{i} = {a_i}")
        print(f"  a_{i+1} = {a_i_plus_1}")
        print(f"  b_{i}^2 = {a_i}^2 + {a_i_plus_1}^2 + 2*{a_i}*{a_i_plus_1}*cos({phi:.4f}) = {b_i_squared:.4f}")
        print(f"  b_{i} = sqrt({b_i_squared:.4f}) = {b_i:.4f}")
        print(f"  h_{i} = ({a_i} * {a_i_plus_1} * sin({phi:.4f})) / (2 * {b_i:.4f}) = {h_i:.4f}")
        print("-" * 20)
        
        if h_i > max_h:
            max_h = h_i
            
    print(f"The values of h_i are distances computed for each vertex corner.")
    print(f"The largest possible Hausdorff distance is the maximum of these values.")
    print(f"\nFinal Result:")
    print(f"The largest possible Hausdorff distance is {max_h}")

    return max_h

if __name__ == '__main__':
    # Example usage: a pentagon with varying edge lengths
    # You can change these values to match your specific polygon
    # a_i = length of edge i
    edge_lengths_example = [10, 12, 8, 15, 11]
    
    # Let's consider a simpler case: a rectangle with sides 6 and 8.
    # The edge lengths would be [6, 8, 6, 8]
    # For a rectangle n=4, phi = pi/2, cos(phi)=0, sin(phi)=1
    # b_i = sqrt(a_i^2 + a_{i+1}^2)
    # h_i = a_i * a_{i+1} / (2 * b_i)
    print("Example: Rectangle with sides 6 and 8")
    rectangle_edges = [6.0, 8.0, 6.0, 8.0]
    calculate_max_hausdorff_distance(rectangle_edges)
    
    print("\n" + "="*50 + "\n")
    
    print("Example: Pentagon with edges [10, 12, 8, 15, 11]")
    pentagon_edges = [10.0, 12.0, 8.0, 15.0, 11.0]
    calculate_max_hausdorff_distance(pentagon_edges)
