def solve():
    """
    This function calculates the coefficients for the number of closed tree-like walks of length 6.
    The formula is:
    c_1 * e + c_2 * k + c_3 * p + c_4 * sum(deg(v) choose 2) + c_5 * sum(deg(v) choose 3)

    The derivation of the coefficients is as follows:
    A closed tree-like walk of length 6 traverses a tree with 3 edges (each edge twice),
    or a tree with 2 edges (one twice, one four times), or a tree with 1 edge (six times).

    1.  c_1: The tree is a single edge (P_2). The walk is u -> v -> u -> v -> u -> v -> u.
        For each edge, there are 2 such walks (one starting at each end). So, c_1 = 2.

    2.  c_2: The subgraph is a triangle (K_3). A triangle is not a tree. By definition, there are no
        tree-like walks whose minimal supporting subgraph is a triangle. So, c_2 = 0.

    3.  c_3: The tree is a path of length 3 (P_4). A detailed enumeration shows there are 6 distinct
        walks that traverse all edges of a P_4 twice. So, c_3 = 6.

    4.  c_4: The tree is a path of length 2 (P_3). One edge is traversed 4 times, the other twice.
        A detailed enumeration shows there are 12 distinct walks for each P_3. So, c_4 = 12.

    5.  c_5: The tree is a star graph with 3 edges (K_{1,3}). All edges are traversed twice.
        A detailed enumeration shows there are 12 distinct walks for each K_{1,3}. So, c_5 = 12.
    """
    c1 = 2
    c2 = 0
    c3 = 6
    c4 = 12
    c5 = 12
    
    print(f"c_1 = {c1}")
    print(f"c_2 = {c2}")
    print(f"c_3 = {c3}")
    print(f"c_4 = {c4}")
    print(f"c_5 = {c5}")

solve()
