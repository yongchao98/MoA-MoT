import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer n-gon approximation B with given edge lengths.

    Args:
        n (int): The number of vertices/edges of the polygon B.
        a (list of float): The lengths of the edges of B, a_1, ..., a_n.
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"The number of edge lengths ({len(a)}) must be equal to n ({n}).")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    print(f"For n = {n}:")
    print(f"Angle phi = 2*pi/n = {phi:.4f} radians")
    print("-" * 30)

    h_i_values = []
    
    # The final formula is D = max_i(H_i/2) where H_i = (a_i * a_{i+1} * sin(phi)) / b_i
    # We first calculate all components H_i and then find the maximum.
    
    print("Calculating components for each vertex i:")
    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic index, a_{n+1} is a_1
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        b_i = math.sqrt(b_i_squared)

        # This is the altitude of the triangle formed by v_{i-1}, v_i, v_{i+1}
        # In the referenced formula, this term is H_i
        H_i = (a_i * a_i_plus_1 * sin_phi) / b_i
        h_i_values.append(H_i)

        print(f"i = {i+1}:")
        print(f"  a_{i} = {a_i}, a_{i+1} = {a_i_plus_1}")
        print(f"  b_{i} = sqrt({a_i**2:.2f} + {a_i_plus_1**2:.2f} + 2*{a_i:.2f}*{a_i_plus_1:.2f}*cos({phi:.2f})) = {b_i:.4f}")
        print(f"  H_{i} = (a_{i} * a_{i+1} * sin(phi)) / b_{i} = {H_i:.4f}")

    # The maximum possible Hausdorff distance is max(H_i) / 2
    max_H = max(h_i_values)
    max_dist = 0.5 * max_H

    print("-" * 30)
    print("The final result is (1/2) * max(H_i)")
    print(f"max(H_i) = {max_H:.4f}")
    print(f"Largest possible Hausdorff distance = 0.5 * {max_H:.4f} = {max_dist:.4f}")
    
    return max_dist

# --- Example Usage ---
# Define the polygon's properties
n_sides = 5
edge_lengths = [3.0, 4.0, 5.0, 6.0, 7.0]

# Calculate and print the result
result = calculate_max_hausdorff_distance(n_sides, edge_lengths)
# The final answer format is specified to be like <<<value>>>
# To adhere to the format, we will output the numerical result from our example case.
print(f"\nFinal Answer in specified format:")
print(f"<<<{result:.4f}>>>")
