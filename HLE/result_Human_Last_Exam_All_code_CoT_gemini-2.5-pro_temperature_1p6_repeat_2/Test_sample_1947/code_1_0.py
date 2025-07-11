def solve():
    """
    This function calculates and prints the coefficients c_1, c_2, c_3, c_4, c_5.

    The number of closed tree-like walks of length 6 in a simple graph X is given by the expression:
    c_1*e + c_2*k + c_3*p + c_4*Sum(deg(v) choose 2) + c_5*Sum(deg(v) choose 3)

    A closed tree-like walk of length 6 is a walk that traverses the 3 edges of a tree subgraph twice (once in each direction).
    The possible tree subgraphs with 3 edges are the path P_4 and the star graph K_{1,3}.

    1. For each subgraph isomorphic to P_4 (a path with 4 vertices), there are 6 such walks.
       Let a path be v1-v2-v3-v4.
       - Starting from endpoints (v1, v4): 2 walks (e.g., v1-v2-v3-v4-v3-v2-v1).
       - Starting from internal vertices (v2, v3): 4 walks (e.g., v2-v1-v2-v3-v4-v3-v2).
       Total per P_4 is 6. The number of P_4 subgraphs is p. Contribution: 6*p.

    2. For each subgraph isomorphic to K_{1,3} (a star graph with one central vertex and 3 leaves), there are 12 such walks.
       Let c be the center and l1, l2, l3 be the leaves.
       - Starting from the center c: 3! = 6 walks (permutations of visiting the leaves, e.g., c-l1-c-l2-c-l3-c).
       - Starting from a leaf (e.g., l1): 2 walks for each of the 3 leaves (e.g., l1-c-l2-c-l3-c-l1). Total: 3 * 2 = 6.
       Total per K_{1,3} is 12. The number of K_{1,3} subgraphs is Sum(deg(v) choose 3). Contribution: 12 * Sum(deg(v) choose 3).

    The total number of such walks is the sum of contributions from all such subgraphs:
    N = 6*p + 12*Sum(deg(v) choose 3).

    Comparing this with the given formula, we get the coefficients:
    c_1 = 0 (no dependency on e)
    c_2 = 0 (no dependency on k)
    c_3 = 6
    c_4 = 0 (no dependency on Sum(deg(v) choose 2))
    c_5 = 12
    """
    
    c1 = 0
    c2 = 0
    c3 = 6
    c4 = 0
    c5 = 12
    
    # Print the coefficients in order
    print(c1)
    print(c2)
    print(c3)
    print(c4)
    print(c5)

solve()