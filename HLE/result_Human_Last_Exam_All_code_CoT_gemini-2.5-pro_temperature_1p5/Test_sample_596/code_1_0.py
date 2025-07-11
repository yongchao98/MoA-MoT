import numpy as np

def solve():
    """
    This function calculates the number of higher dimensional rooted forests (F,R)
    of the standard triangulation of the Möbius band that fail to have the forest F
    simplicially collapse onto the root R.

    The method is based on formulas from algebraic combinatorics:
    - Total number of rooted 2-forests (F,R) is det(2I + d2^T * d2).
    - Number of rooted 2-forests (F,R) where F collapses to R is det(I + d2^T * d2).
    - The difference gives the number of non-collapsing pairs.

    We use a minimal triangulation of the Möbius band with:
    - 5 vertices: {1, 2, 3, 4, 5}
    - 10 edges (all pairs of vertices, the graph K5)
    - 5 faces: {{1,2,5}, {1,3,4}, {2,3,5}, {1,2,4}, {3,4,5}}
    """

    # Define the oriented edges and faces of the triangulation.
    # Edges are (i, j) with i < j.
    edges = [(1,2), (1,3), (1,4), (1,5), (2,3), (2,4), (2,5), (3,4), (3,5), (4,5)]
    edge_to_idx = {edge: i for i, edge in enumerate(edges)}
    num_edges = len(edges)

    # Faces are sets of 3 vertices. We'll orient them as (i,j,k) with i<j<k.
    faces_set = [{1,2,5}, {1,3,4}, {2,3,5}, {1,2,4}, {3,4,5}]
    faces = [tuple(sorted(list(f))) for f in faces_set]
    face_to_idx = {face: i for i, face in enumerate(faces)}
    num_faces = len(faces)
    
    # Initialize the boundary matrix d2 (f1 x f2)
    d2 = np.zeros((num_edges, num_faces), dtype=int)

    # Populate the boundary matrix d2
    # For a face (i,j,k) with i<j<k, the oriented boundary is (i,j)-(i,k)+(j,k)
    for face_idx, face in enumerate(faces):
        i, j, k = face
        
        # edge (i,j)
        e1_idx = edge_to_idx[(i,j)]
        d2[e1_idx, face_idx] += 1
        
        # edge (i,k)
        e2_idx = edge_to_idx[(i,k)]
        d2[e2_idx, face_idx] -= 1
        
        # edge (j,k)
        e3_idx = edge_to_idx[(j,k)]
        d2[e3_idx, face_idx] += 1

    # Compute M = d2^T * d2
    M = d2.T @ d2
    
    # The identity matrix of the same size as M
    I = np.identity(num_faces)
    
    # Calculate the two determinants. These will be integers.
    det_collapsing = np.linalg.det(I + M)
    det_total = np.linalg.det(2*I + M)
    
    # The result must be an integer, so we round the floating point results.
    num_collapsing = int(round(det_collapsing))
    num_total = int(round(det_total))
    
    # The number of non-collapsing pairs is the difference.
    num_non_collapsing = num_total - num_collapsing

    print(f"The total number of rooted 2-forests is {num_total}.")
    print(f"The number of rooted 2-forests that simplicially collapse is {num_collapsing}.")
    print("The number of forests that fail to collapse is the difference:")
    print(f"{num_total} - {num_collapsing} = {num_non_collapsing}")

solve()
<<<175>>>