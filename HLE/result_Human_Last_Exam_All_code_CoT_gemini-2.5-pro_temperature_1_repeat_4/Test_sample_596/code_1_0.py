import numpy as np
from sympy import Matrix

def compute_homology_torsion():
    """
    Computes the order of the torsion subgroup of the first homology group
    of a standard triangulation of the Möbius band.
    """
    # Standard combinatorial triangulation of the Mobius band (mb_7_10)
    # 7 vertices, 17 edges, 10 faces
    vertices = list(range(1, 8))
    
    # Edges are defined as sorted tuples of vertices
    edges = sorted(list(set(
        tuple(sorted((f[i], f[j])))
        for f in [
            (1,2,4), (1,2,5), (1,3,6), (1,3,7), (1,4,5), (1,6,7), 
            (2,3,4), (2,3,7), (2,5,6), (2,6,7)
        ]
        for i in range(3) for j in range(i + 1, 3)
    )))
    
    # Faces are defined as sorted tuples of vertices
    faces = [
        (1,2,4), (1,2,5), (1,3,6), (1,3,7), (1,4,5), (1,6,7), 
        (2,3,4), (2,3,7), (2,5,6), (2,6,7)
    ]

    num_vertices = len(vertices)
    num_edges = len(edges)
    num_faces = len(faces)

    vert_map = {v: i for i, v in enumerate(vertices)}
    edge_map = {e: i for i, e in enumerate(edges)}

    # Build boundary matrix d1: C1 -> C0 (edges to vertices)
    d1 = np.zeros((num_vertices, num_edges), dtype=int)
    for j, edge in enumerate(edges):
        v1, v2 = edge
        d1[vert_map[v2], j] = 1
        d1[vert_map[v1], j] = -1

    # Build boundary matrix d2: C2 -> C1 (faces to edges)
    d2 = np.zeros((num_edges, num_faces), dtype=int)
    for j, face in enumerate(faces):
        v1, v2, v3 = sorted(face)
        # Orient face boundary as (v1,v2) + (v2,v3) - (v1,v3)
        d2[edge_map[tuple(sorted((v1, v2)))], j] = 1 if v1 < v2 else -1
        d2[edge_map[tuple(sorted((v2, v3)))], j] = 1 if v2 < v3 else -1
        d2[edge_map[tuple(sorted((v1, v3)))], j] = -1 if v1 < v3 else 1
    
    # Use sympy for Smith Normal Form calculation
    d1_matrix = Matrix(d1)
    d2_matrix = Matrix(d2)

    # Algorithm to find torsion of H1
    # 1. Compute SNF of d1: U1 * d1 * V1 = D1
    D1, U1, _ = d1_matrix.smith_normal_form(domain='Z')
    
    # 2. Transform d2 using U1: d2_hat = U1 * d2
    d2_hat = U1 * d2_matrix
    
    # 3. Find torsion coefficients
    torsion_coeffs = []
    # Diagonal entries of D1 are the invariant factors of the domain of d1
    d1_invariants = D1.diagonal()
    
    # Rank of d1
    rank_d1 = sum(1 for d in d1_invariants if d != 0)

    # Torsion appears where d1 has non-unit invariant factors
    # and the corresponding cycles are not boundaries from d2.
    # This corresponds to rows in d2_hat being zero for a non-unit invariant factor in D1.
    for i in range(rank_d1):
        d = d1_invariants[i]
        if d > 1:
            # Check if the i-th row of d2_hat is zero
            is_boundary = any(d2_hat[i, j] != 0 for j in range(d2_hat.cols))
            if not is_boundary:
                torsion_coeffs.append(d)

    if not torsion_coeffs:
        print("The first homology group has no torsion.")
        print("The number of such forests is the order of the trivial group.")
        print("1")
        return 1
    else:
        # Calculate the order of the torsion subgroup
        order = 1
        for coeff in torsion_coeffs:
            order *= coeff
        
        equation_parts = [str(c) for c in torsion_coeffs]
        print(f"The torsion coefficients are: {torsion_coeffs}")
        print("The number of such forests is the product of these coefficients.")
        print(f"{' * '.join(equation_parts)} = {order}")
        return order

if __name__ == '__main__':
    compute_homology_torsion()
