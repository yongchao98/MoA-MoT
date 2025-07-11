def solve_coefficients():
    """
    This function determines the coefficients c_1, c_2, c_3, c_4, c_5 for the number
    of closed tree-like walks of length 6 in a simple graph X.

    The reasoning is based on combinatorial counting of walk patterns on small tree subgraphs.
    A closed tree-like walk of length 6 must have an underlying tree structure with 1, 2, or 3 edges.

    Let's determine each coefficient by counting the number of possible walks for each base subgraph.
    """

    # c1: Underlying tree is P2 (one edge).
    # The walk must traverse the edge 6 times (3 times out, 3 times back).
    # For an edge (u,v), walks can start at u or v.
    # Walk starting at u: (u,v,u,v,u,v,u) -> 1 walk
    # Walk starting at v: (v,u,v,u,v,u,v) -> 1 walk
    # Total walks per edge = 2.
    c1 = 2

    # c2: For subgraphs isomorphic to K3 (triangle).
    # A "tree-like" walk is one where the set of edges traversed forms a tree.
    # A K3 is a cycle, not a tree. Any tree-like walk on a graph with K3's must use a
    # subset of the K3's edges that forms a tree (i.e., P2 or P3).
    # These cases are already counted by other terms (c1*e and c4*sum(deg choose 2)).
    # Thus, there's no unique contribution from K3 subgraphs.
    c2 = 0

    # c3: Underlying tree is P4 (path of 3 edges).
    # The walk traverses each of the 3 edges twice.
    # For a path v1-v2-v3-v4:
    # Walks starting at endpoints (v1, v4): Go to the other end and back. 1 walk each. (2 total)
    # Walks starting at interior vertices (v2, v3): Perform a round trip on each side.
    # The order can be varied. e.g., from v2: (v2->v1->v2) then (v2->v3->v4->v3->v2) or vice versa. 2 walks each. (4 total)
    # Total walks per P4 = 2 + 4 = 6.
    c3 = 6

    # c4: Underlying tree is P3 (path of 2 edges).
    # One edge is traversed 4 times, the other 2 times. Let path be u-v-w.
    # Case 1: edge (u,v) 4x, (v,w) 2x.
    # Starts at v (center): 3 permutations of loops. 3 walks.
    # Starts at u (leaf): 2 ways to order the inner loops. 2 walks.
    # Starts at w (other leaf): 1 walk (inner part must use (u,v) 4x).
    # Total = 3 + 2 + 1 = 6 walks.
    # Case 2: edge (u,v) 2x, (v,w) 4x. Symmetric, so 6 walks.
    # Total walks per P3 = 6 + 6 = 12.
    c4 = 12

    # c5: Underlying tree is K1,3 (star graph with 3 edges).
    # Each of the 3 edges is traversed twice. Center v, leaves u1, u2, u3.
    # Starts at v (center): Permutations of the 3 loops v-ui-v. 3! = 6 walks.
    # Starts at a leaf (e.g. u1): A round trip u1-v-u1 plus an inner walk from v using u2,u3.
    # Inner walk has 2 permutations. 2 walks per leaf. 3*2=6 walks.
    # Total walks per K1,3 = 6 + 6 = 12.
    c5 = 12

    print(f"{c1},{c2},{c3},{c4},{c5}")

solve_coefficients()