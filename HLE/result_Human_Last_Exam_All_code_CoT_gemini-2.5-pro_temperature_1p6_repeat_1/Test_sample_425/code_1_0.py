import numpy as np

def demonstrate_projections():
    """
    This function demonstrates how planar projections of a 3D object
    with A4 rotational symmetry can have rotational symmetries of order 3 and 4.
    """
    print("--- Demonstration for Possible Projection Symmetries ---")

    # Case i): Order 3
    # We construct a tetrahedron with a 3-fold axis along the z-axis.
    # Its projection onto the xy-plane will have 3-fold symmetry.
    # The vertices are placed symmetrically around the z-axis.
    r_3 = np.sqrt(8/9)
    v1_3 = np.array([0, 0, 1])
    v2_3 = np.array([r_3, 0, -1/3.0])
    v3_3 = np.array([r_3 * np.cos(2*np.pi/3), r_3 * np.sin(2*np.pi/3), -1/3.0])
    v4_3 = np.array([r_3 * np.cos(4*np.pi/3), r_3 * np.sin(4*np.pi/3), -1/3.0])
    A_for_C3 = [v1_3, v2_3, v3_3, v4_3]
    B_proj_C3 = [p[:2] for p in A_for_C3]
    
    print("\ni) Possibility of Order 3:")
    print("A is a set of 4 points forming a regular tetrahedron.")
    print("Projection B is onto the xy-plane along a 3-fold symmetry axis.")
    print("The projected points are:")
    for p in B_proj_C3:
        print(f"  ({p[0]:.4f}, {p[1]:.4f})")
    print("These points (one at the center, three forming an equilateral triangle) have C3 rotational symmetry (order 3).")

    # Case ii): Order 4
    # We construct a tetrahedron whose vertices are selected from the vertices of a cube.
    # This specific orientation, when projected along a coordinate axis, results in a square.
    # This is an example of "accidental symmetry".
    v1_4 = np.array([1, 1, 1])
    v2_4 = np.array([1, -1, -1])
    v3_4 = np.array([-1, 1, -1])
    v4_4 = np.array([-1, -1, 1])
    A_for_C4 = [v1_4, v2_4, v3_4, v4_4]
    B_proj_C4 = [p[:2] for p in A_for_C4]

    print("\nii) Possibility of Order 4:")
    print("A is a set of 4 points, which also form a regular tetrahedron.")
    print("Its rotational symmetry group is A4.")
    print("Projection B is onto the xy-plane along the z-axis.")
    print("The projected points are:")
    for p in B_proj_C4:
        print(f"  ({p[0]:.4f}, {p[1]:.4f})")
    print("These four points are the vertices of a square, which has C4 rotational symmetry (order 4).")
    
    print("\niii) and iv) Order 6 and Infinity are not possible:")
    print("These would require the 3D object A to have symmetries (e.g., 6-fold or continuous rotation) that are not part of the A4 group.")

demonstrate_projections()