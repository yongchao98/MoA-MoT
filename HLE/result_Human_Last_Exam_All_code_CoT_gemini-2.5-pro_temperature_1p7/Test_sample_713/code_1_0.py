import numpy as np
from scipy.spatial import ConvexHull

def solve_polytope_volume():
    """
    Demonstrates that for a given d-simplex T with volume V,
    the volume of the associated polytope P is d! * V.

    This function uses a specific 3D simplex as an example.
    It requires the 'pycddlib' library for vertex enumeration
    from hyperplane representations. You can install it using:
    pip install pycddlib
    """
    try:
        import cdd
    except ImportError:
        print("This script requires the 'pycddlib' library.")
        print("Please install it by running: pip install pycddlib")
        return

    # 1. Define a sample 3D simplex T
    d = 3
    # Vertices of the simplex T. v0 is at the origin.
    v0 = np.array([0, 0, 0])
    v1 = np.array([2, 0, 0])
    v2 = np.array([1, 2, 0])
    v3 = np.array([1, 1, 2])
    vertices = [v0, v1, v2, v3]

    # 2. Calculate the volume of the simplex T
    # For a simplex with one vertex at the origin, V = (1/d!) * |det(v1, v2, ...)|
    M = np.array([v1, v2, v3]).T
    simplex_volume_V = np.abs(np.linalg.det(M)) / np.math.factorial(d)
    
    # 3. Formulate the linear inequalities for the polytope P
    # The polytope P is defined by strips for each edge e_ij = v_j - v_i.
    # Each strip is c_1 <= n . x <= c_2, which gives two inequalities:
    # n . x <= c_2 and -n . x <= -c_1.
    # The format for cddlib is Ax + b >= 0, which is -A'x + b' >= 0
    # where Ax <= b becomes A'x >= b', and our b vector is on the right side of the matrix.
    # Matrix H = [b, -A]
    inequalities = []
    edges = [
        (v1 - v0, v0, v1), (v2 - v0, v0, v2), (v3 - v0, v0, v3),
        (v2 - v1, v1, v2), (v3 - v1, v1, v3), (v3 - v2, v2, v3)
    ]

    for edge_vec, start_v, end_v in edges:
        c1 = np.dot(edge_vec, start_v)
        c2 = np.dot(edge_vec, end_v)
        
        # Ensure c1 <= c2
        if c1 > c2:
            c1, c2 = c2, c1
        
        # n . x <= c2  => -n . x + c2 >= 0
        inequalities.append([c2] + list(-edge_vec))
        # n . x >= c1  => n . x - c1 >= 0
        inequalities.append([-c1] + list(edge_vec))

    # 4. Compute vertices of P using pycddlib
    # Create the cdd matrix representation for H-representation
    h_rep = cdd.Matrix(inequalities, number_type='float')
    h_rep.rep_type = cdd.RepType.INEQUALITY
    
    # Create the Polyhedron object
    poly = cdd.Polyhedron(h_rep)
    
    # Get the V-representation (vertices)
    v_rep = poly.get_generators()
    polytope_vertices = np.array(v_rep)
    
    # The output includes vertices (linearity=0) and rays (linearity=1).
    # For a bounded polytope, there are no rays. We filter for vertices.
    if polytope_vertices.shape[1] > d + 1: # check if rays are present
       polytope_vertices = polytope_vertices[polytope_vertices[:,0] == 1][:, 1:]
    else: # old pycddlib version
       polytope_vertices = polytope_vertices[:, 1:]

    # 5. Calculate the volume of P
    hull = ConvexHull(polytope_vertices)
    polytope_volume_P = hull.volume

    # 6. Print the results
    d_factorial_V = np.math.factorial(d) * simplex_volume_V
    
    print(f"Dimension d = {d}")
    print(f"Simplex T volume V = {simplex_volume_V:.4f}")
    print(f"Polytope P volume (calculated) = {polytope_volume_P:.4f}")
    print(f"Theoretical P volume ({d}! * V) = {d_factorial_V:.4f}")
    print("\nVerification:")
    print(f"Is Vol(P) ≈ {d}! * V? {'Yes' if np.isclose(polytope_volume_P, d_factorial_V) else 'No'}")

    # Outputting the final equation as requested
    print("\nThe relationship is expressed as:")
    print(f"Vol(P) = {d}! x V")
    print(f"{polytope_volume_P:.4f} = {np.math.factorial(d)} x {simplex_volume_V:.4f}")

if __name__ == '__main__':
    solve_polytope_volume()
