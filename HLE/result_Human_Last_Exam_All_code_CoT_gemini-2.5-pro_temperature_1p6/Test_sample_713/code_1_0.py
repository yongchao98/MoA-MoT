import numpy as np
from scipy.spatial import ConvexHull
# The pypoman library is required for converting hyperplane representation to vertex representation.
# You can install it using: pip install pypoman
import pypoman

def get_simplex_volume(vertices, d):
    """Calculates the volume of a simplex given its vertices."""
    # We take d vectors from one vertex (e.g., v0) to the others.
    v0 = vertices[0]
    vectors = vertices[1:] - v0
    # The volume is 1/d! times the volume of the parallelepiped
    # spanned by these vectors.
    M = np.array(vectors).T
    volume = np.abs(np.linalg.det(M)) / np.math.factorial(d)
    return volume

def get_polytope_volume_from_H_representation(A, b):
    """
    Calculates the volume of a polytope given its H-representation (Ax <= b).
    It first computes the vertices of the polytope.
    """
    try:
        # Use pypoman to find the vertices of the polytope.
        vertices = pypoman.compute_polytope_vertices(A, b)
    except Exception as e:
        print(f"Could not compute vertices: {e}")
        return 0

    if not vertices:
        print("Polytope is empty or unbounded.")
        return 0
    
    # Ensure we have enough points to form a hull in d-dimensions.
    if len(vertices) < d + 1:
        return 0

    try:
        # The volume of the polytope is the volume of the convex hull of its vertices.
        hull = ConvexHull(np.array(vertices))
        return hull.volume
    except Exception as e:
        print(f"Could not compute convex hull: {e}")
        return 0


def run_verification(name, vertices):
    """
    Runs the verification for a given simplex.
    """
    print(f"--- {name} ---")
    d = vertices.shape[1]
    num_vertices = d + 1
    num_edges = num_vertices * (num_vertices - 1) // 2

    # 1. Calculate the volume of the simplex T
    V_T = get_simplex_volume(vertices, d)
    if V_T < 1e-9:
        print("The simplex is degenerate (volume is near zero).")
        return

    # 2. Construct the H-representation (Ax <= b) of the polytope P
    A = np.zeros((num_edges * 2, d))
    b = np.zeros(num_edges * 2)
    
    row_idx = 0
    for i in range(num_vertices):
        for j in range(i + 1, num_vertices):
            v_i = vertices[i]
            v_j = vertices[j]
            edge_vec = v_j - v_i
            
            # The two hyperplanes for the edge (v_i, v_j) define a slab.
            # edge_vec . x <= edge_vec . v_j
            # -edge_vec . x <= -edge_vec . v_i
            A[row_idx] = edge_vec
            b[row_idx] = np.dot(edge_vec, v_j)
            row_idx += 1
            
            A[row_idx] = -edge_vec
            b[row_idx] = -np.dot(edge_vec, v_i)
            row_idx += 1

    # 3. Calculate the volume of the polytope P
    V_P = get_polytope_volume_from_H_representation(A, b)

    # 4. Print results and the final equation
    print(f"Dimension d = {d}")
    print(f"Volume of simplex T (V) = {V_T:.6f}")
    print(f"Volume of polytope P = {V_P:.6f}")
    
    ratio = V_P / V_T
    d_factorial = float(np.math.factorial(d))
    
    print("\nFinal Equation:")
    print(f"Vol(P) / Vol(T) = {V_P:.6f} / {V_T:.6f} = {ratio:.6f}")
    print(f"This should be equal to d! = {d}! = {d_factorial}")
    print("-" * (len(name) + 8) + "\n")


if __name__ == '__main__':
    # --- Case d=2 ---
    # Example 1: Right-angled triangle in R^2
    vertices_T1 = np.array([[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]])
    run_verification("Case 1: Right-angled triangle in R^2", vertices_T1)

    # Example 2: Equilateral-like triangle in R^2
    l = 4.0
    h = l * np.sqrt(3)/2
    vertices_T2 = np.array([[0,0], [l, 0], [l/2, h]])
    run_verification("Case 2: Equilateral triangle in R^2", vertices_T2)

    # Example 3: A generic 'thin' triangle in R^2
    vertices_T3 = np.array([[0, 0], [10, 1], [5, 8]])
    run_verification("Case 3: A generic triangle in R^2", vertices_T3)
    
    # --- Case d=3 ---
    # Example 4: An orthogonal simplex (corner of a box) in R^3
    vertices_T4 = np.array([[0.,0.,0.], [2.,0.,0.], [0.,3.,0.], [0.,0.,4.]])
    run_verification("Case 4: Orthogonal simplex in R^3", vertices_T4)
    
    # Example 5: A random simplex in R^3
    # Use a fixed seed for reproducibility
    np.random.seed(0)
    vertices_T5 = np.random.rand(4, 3) * 5
    run_verification("Case 5: A random simplex in R^3", vertices_T5)
