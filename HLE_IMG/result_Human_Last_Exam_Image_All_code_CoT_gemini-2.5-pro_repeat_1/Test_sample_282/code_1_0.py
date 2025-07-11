import numpy as np

def solve():
    """
    This script evaluates statement D by calculating and comparing the number of walks
    in two different scenarios described in the statement.
    """
    num_vertices = 9
    
    # 1. Define the adjacency matrix for the directed graph from the image.
    adj_list = {
        0: [3], 1: [8], 2: [4], 3: [2, 6], 4: [1, 2, 5], 
        5: [8], 6: [1, 5], 7: [3], 8: [7]
    }
    A_dir = np.zeros((num_vertices, num_vertices), dtype=int)
    for u, neighbors in adj_list.items():
        for v in neighbors:
            A_dir[u, v] = 1

    # 2. Calculate the number of walks from vertex 3 to 7 with length k <= 10 in the digraph.
    total_walks_digraph = 0
    Ak_dir = np.copy(A_dir)
    for k in range(1, 11):
        if k > 1:
            Ak_dir = np.dot(Ak_dir, A_dir)
        total_walks_digraph += Ak_dir[3, 7]

    # 3. Create the adjacency matrix for the undirected graph.
    # An edge (u,v) exists if there is an arc u->v or v->u in the digraph.
    A_undir = np.logical_or(A_dir, A_dir.T).astype(int)

    # 4. Calculate the number of walks from vertex 6 to 5 with length k <= 3 in the undirected graph.
    total_walks_undir = 0
    Ak_undir = np.copy(A_undir)
    for k in range(1, 4):
        if k > 1:
            Ak_undir = np.dot(Ak_undir, A_undir)
        total_walks_undir += Ak_undir[6, 5]

    # 5. Print the results and the conclusion.
    print("Evaluating Statement D:")
    print(f"The number of walks from vertex 3 to vertex 7 with length k <= 10 is: {total_walks_digraph}")
    print(f"The number of walks from vertex 6 to vertex 5 with length k <= 3 is: {total_walks_undir}")
    
    # The final equation from the statement
    print("\nFinal equation check:")
    print(f"{total_walks_digraph} == {total_walks_undir}")

    if total_walks_digraph == total_walks_undir:
        print("\nSince the two quantities are equal, statement D is correct.")
    else:
        print("\nSince the two quantities are not equal, statement D is incorrect.")

solve()