import numpy as np
from sympy import Matrix

def solve_mobius_forest_problem():
    """
    This function computes the number of non-collapsing higher dimensional rooted forests
    on the standard triangulation of the Möbius band, interpreted as the size of the
    torsion subgroup of its critical group.
    """
    
    # We use the M_3 triangulation of the Möbius band. It has 6 vertices, 12 edges, and 6 faces.
    # We define an orientation and construct the boundary matrix d2, mapping C_2 (faces) to C_1 (edges).
    # The matrix below represents the boundary map from the 6 faces to the 12 edges.
    # Each column corresponds to a face, and each row to an edge.
    # The choice of basis and orientation is consistent with standard literature.
    
    # Edges (rows): e1..e12. Faces (columns): f1..f6.
    d2_np = np.array([
        [ 0,  1,  0,  0,  0,  0],  # e1
        [-1,  0,  0,  0,  0,  0],  # e2
        [ 0,  0,  1,  0,  0,  0],  # e3
        [-1,  0,  0,  1,  0,  0],  # e4
        [ 0, -1,  0, -1,  0,  0],  # e5
        [ 0,  1,  0,  0,  1,  0],  # e6
        [ 0,  0, -1,  0, -1,  0],  # e7
        [ 1,  0,  0,  0,  0,  1],  # e8
        [ 0,  0,  1,  0,  0, -1],  # e9
        [ 0,  0,  0,  1,  0,  0],  # e10
        [ 0,  0,  0,  0,  0,  1],  # e11
        [ 0,  0,  0,  0,  1,  0]   # e12
    ])

    # The torsion part of the critical group is determined by the matrix L = d2^T * d2.
    L_np = d2_np.T @ d2_np
    
    # We use sympy to compute the Smith Normal Form of L over the integers.
    L_sympy = Matrix(L_np)
    snf_matrix, _, _ = L_sympy.smith_normal_form()
    
    # The diagonal elements of the SNF matrix are the invariant factors.
    invariant_factors = snf_matrix.diagonal()
    
    # The size of the torsion group is the product of the invariant factors greater than 1.
    torsion_factors = [int(f) for f in invariant_factors if f > 1]
    
    result = 1
    if not torsion_factors:
        result = 1 # Corresponds to trivial torsion group
    else:
        for factor in torsion_factors:
            result *= factor

    print("The number of higher dimensional rooted forests on the standard triangulation of the Möbius band that fail to simplicially collapse onto the root is interpreted as the size of the torsion part of its critical group.")
    print("This is calculated from the Smith Normal Form of the matrix L = d2^T * d2.")
    print(f"The invariant factors of L are: {[int(f) for f in invariant_factors]}")
    print(f"The non-trivial invariant factors (those > 1) are: {torsion_factors}")

    if len(torsion_factors) > 1:
        equation_str = " * ".join(map(str, torsion_factors))
        print(f"The final result is the product of these factors: {equation_str} = {result}")
    else:
        print(f"The final result is: {result}")

solve_mobius_forest_problem()