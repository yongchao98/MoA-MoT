import sys

def solve():
    """
    This function solves the problem about the minimal number of edges
    to add to G' to make it 2-edge-connected.

    The reasoning is as follows:
    1. The number of edges to add to a graph to make it 2-edge-connected is determined
       by its structure, specifically the number of "leaf blocks" (k) in its
       block-cut forest. The formula is ceil(k/2).
    2. The question asks for the "minimal number", which implies we are looking for the
       minimum possible value of this quantity across all valid graphs G.
    3. If we can construct a graph G satisfying all the conditions such that the resulting
       graph G' is already 2-edge-connected, then for this G', k=0. The number of
       edges to add would be ceil(0/2) = 0.
    4. Since the number of added edges cannot be negative, 0 is the minimum possible value.
    5. We need to demonstrate that such a G can exist. Let's take d=2 (an even integer).
       The degrees of v1, v2, v3 are 2, 3, and 3. Let G' be a cycle graph C_n (e.g., n=3),
       which is 2-edge-connected.
       We can form a graph G with edge connectivity 2 by connecting all edges from v1, v2, v3
       to the vertices of G'. For instance, connect all of them to a single vertex of C_n.
       It can be verified that the resulting graph G has an edge connectivity of 2.
    6. Since it's possible to have a G' that is already 2-edge-connected, the minimal
       number of edges one must add is 0.

    The variables in the problem are d, but the result is a constant value independent of d.
    The final equation is simply the number 0.
    """
    
    # The minimal number of edges is 0.
    minimal_edges = 0
    
    # We need to output the numbers in the final equation.
    # The final equation is just the answer itself.
    print(minimal_edges)

solve()
