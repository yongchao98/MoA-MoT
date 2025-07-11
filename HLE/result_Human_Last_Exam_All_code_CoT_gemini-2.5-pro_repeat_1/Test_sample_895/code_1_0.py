import math

def solve_hausdorff_distance():
    """
    This script calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B, given the edge lengths of B.
    """
    
    # --- Problem Parameters ---
    # You can modify these values to solve for a different polygon.
    # n: number of sides of the polygon B.
    # a: list of edge lengths of B. The indexing follows the problem:
    #    vertex v_i is the intersection of the lines defining edges E_i and E_{i+1}.
    #    The length of edge E_i is a_i.
    
    # Example 1: An equilateral-like triangle (n=3)
    n = 3
    a = [10.0, 3.0, 10.0]

    # Example 2: A regular pentagon (n=5)
    # n = 5
    # a = [1.0, 1.0, 1.0, 1.0, 1.0]

    print("--- Problem Setup ---")
    print(f"Number of sides n = {n}")
    print(f"Edge lengths a = {a}\n")

    if n < 3:
        print("Error: A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print("Error: The number of edge lengths must match n.")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_dist = -1.0
    best_case_info = {}

    # Iterate through each vertex v_i of the polygon B
    for i in range(n):
        a_i = a[i]
        a_i_plus_1 = a[(i + 1) % n]
        
        current_dist = 0
        case_type = ""
        params = {}

        # Check for obtuse angle cases. An obtuse angle in the triangle
        # Delta(v_i, v_{i-1}, v_{i+1}) at v_{i+1} or v_{i-1} means the
        # closest point to v_i on the segment [v_{i-1}, v_{i+1}] is an endpoint.
        
        # Case 1: Angle at v_{i+1} is obtuse/right. Distance is d(v_i, v_{i+1}) = a_{i+1}
        if a_i_plus_1 + a_i * cos_phi <= 1e-9: # Use tolerance for float comparison
            current_dist = a_i_plus_1
            case_type = "obtuse"
            params = {'a_i': a_i, 'a_i+1': a_i_plus_1, 'dist': a_i_plus_1}
        # Case 2: Angle at v_{i-1} is obtuse/right. Distance is d(v_i, v_{i-1}) = a_i
        elif a_i + a_i_plus_1 * cos_phi <= 1e-9:
            current_dist = a_i
            case_type = "obtuse"
            params = {'a_i': a_i, 'a_i+1': a_i_plus_1, 'dist': a_i}
        # Case 3: Both base angles are acute. Distance is the altitude h_i.
        else:
            case_type = "acute"
            b_i_sq = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
            b_i = math.sqrt(max(0, b_i_sq))
            
            if b_i < 1e-9:
                h_i = 0
            else:
                h_i = (a_i * a_i_plus_1 * sin_phi) / b_i
            current_dist = h_i
            params = {'a_i': a_i, 'a_i+1': a_i_plus_1, 'phi': phi, 'b_i': b_i, 'h_i': h_i}

        if current_dist > max_dist:
            max_dist = current_dist
            best_case_info = {
                "max_dist": max_dist,
                "index": i,
                "case": case_type,
                "params": params
            }

    # --- Final Result ---
    i = best_case_info['index']
    case = best_case_info['case']
    params = best_case_info['params']
    a_i_val = params.get('a_i')
    a_i_plus_1_val = params.get('a_i+1')
    
    print("--- Calculation of the Maximum Distance ---")
    print(f"The largest possible Hausdorff distance is: {best_case_info['max_dist']:.6f}")
    print(f"\nThis maximum is found by considering the vertex v_{i}, where i={i}.")
    print(f"This vertex corresponds to an '{case}' case.")

    if case == 'acute':
        phi_val = params.get('phi')
        b_i_val = params.get('b_i')
        h_i_val = params.get('h_i')
        print(f"The distance is the altitude h_{i} of the triangle formed by v_{i} and its neighbors.")
        print(f"The calculation is: h_{i} = (a_{i} * a_{{{i}+1}} * sin(phi)) / b_{i}")
        print( "where:")
        print(f"  a_{i} = {a_i_val}")
        print(f"  a_{{{i}+1}} = {a_i_plus_1_val}")
        print(f"  phi = 2*pi/{n} = {phi_val:.6f} rad")
        print(f"  b_{i} = sqrt({a_i_val}^2 + {a_i_plus_1_val}^2 + 2*{a_i_val}*{a_i_plus_1_val}*cos({phi_val:.6f})) = {b_i_val:.6f}")
        print(f"Result: h_{i} = ({a_i_val} * {a_i_plus_1_val} * {math.sin(phi_val):.6f}) / {b_i_val:.6f} = {h_i_val:.6f}")
    
    elif case == 'obtuse':
        dist_val = params.get('dist')
        print("The distance is the length of one of the adjacent edges, as the angle in the triangle is obtuse.")
        print( "The values for this vertex are:")
        print(f"  a_{i} = {a_i_val}")
        print(f"  a_{{{i}+1}} = {a_i_plus_1_val}")
        print(f"Resulting Distance = {dist_val:.6f}")

solve_hausdorff_distance()