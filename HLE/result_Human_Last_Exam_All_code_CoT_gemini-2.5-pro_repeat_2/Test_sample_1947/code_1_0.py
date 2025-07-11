def solve():
    """
    This function determines the coefficients for the number of closed tree-like walks of length 6.

    Let N6_tree be the number of closed tree-like walks of length 6.
    The problem states that N6_tree can be written as:
    c1*e + c2*k + c3*p + c4*Sum(deg(v) choose 2) + c5*Sum(deg(v) choose 3)

    A closed walk of length 6 is "tree-like" if the set of its distinct edges forms a tree.
    This implies the walk is formed by traversing edges of a tree and backtracking.
    The total length of such a walk on a tree T is 2 * (sum of edges traversed).
    If the walk has length 6, the sum of edges traversed (counting multiplicities) is 3.
    Let T be the underlying tree of the walk. Let k_i be the number of times edge e_i of T is traversed back and forth.
    The length of the walk is sum(2*k_i) = 6, so sum(k_i) = 3.

    We classify walks by the number of edges in their underlying tree.
    """

    # Case 1: The underlying tree has 1 edge (P_2).
    # sum(k_i) = k_1 = 3. One edge is traversed back and forth 3 times.
    # For an edge (u, v), the walks are u->v->u->v->u->v->u and v->u->v->u->v->u->v.
    # There are 2 walks per edge. This corresponds to the term 'e'.
    c1 = 2

    # Case 2: The underlying tree is a K_3 (triangle).
    # A K_3 is not a tree. Any walk whose minimal underlying graph is a K_3 is not tree-like.
    # Tree-like walks on graphs containing triangles are counted by other terms based on their
    # actual tree structure (e.g., P_3). So, the coefficient for k is 0.
    c2 = 0

    # Case 3: The underlying tree has 3 edges.
    # sum(k_i) = k_1+k_2+k_3 = 3. Since all edges are used, k_1=k_2=k_3=1.
    # The tree can be a P_4 or a K_{1,3}.

    # Subcase 3.1: The tree is a P_4 (path with 4 vertices, 3 edges).
    # This corresponds to the term 'p'.
    # Let the path be a-b-c-d. Each edge is traversed once back and forth.
    # Walks starting at endpoints (a, d): 1 each (e.g., a->b->c->d->c->b->a). Total 2.
    # Walks starting at inner vertices (b, c): 2 each (e.g., b->a->b->c->d->c->b). Total 4.
    # Total walks per P_4 is 2 + 4 = 6.
    c3 = 6

    # Subcase 3.2: The tree is a K_{1,3} (star graph with 3 edges).
    # This corresponds to Sum(deg(v) choose 3), which counts K_{1,3} subgraphs.
    # Let the center be c, leaves l1, l2, l3.
    # Walk starts at center c: 3! = 6 permutations of excursions c->li->c.
    # Walk starts at a leaf li: The walk is li->c...c->li. The middle part consists of
    # visiting the other 2 leaves from c, which gives 2! = 2 permutations.
    # There are 3 leaves, so 3 * 2 = 6 walks.
    # Total walks per K_{1,3} is 6 + 6 = 12.
    c5 = 12

    # Case 4: The underlying tree has 2 edges (P_3).
    # This corresponds to Sum(deg(v) choose 2), which counts P_3 subgraphs.
    # Let the path be u-v-w. Edges are e1=(u,v), e2=(v,w).
    # sum(k_i) = k_1 + k_2 = 3. Possible (k_1, k_2) are (1,2) and (2,1).
    # Case (k_1,k_2) = (2,1): 2 round trips on e1, 1 on e2.
    #   - Start u: 2 walks.
    #   - Start w: 1 walk.
    #   - Start v: 3!/2! = 3 walks.
    #   Total = 2+1+3 = 6 walks.
    # Case (k_1,k_2) = (1,2): Symmetric, 6 walks.
    # Total walks per P_3 is 6 + 6 = 12.
    c4 = 12

    # The coefficients c1, c2, c3, c4, c5 in order.
    coefficients = [c1, c2, c3, c4, c5]
    
    # The final equation is:
    # 2 * e + 0 * k + 6 * p + 12 * Sum(deg(v) choose 2) + 12 * Sum(deg(v) choose 3)
    # The problem asks for the coefficients c_1, c_2, c_3, c_4, c_5 in this order.
    # I am printing them with some context.
    
    print(f"The number of closed tree-like walks of length 6 is:")
    print(f"{c1} * e + {c2} * k + {c3} * p + {c4} * Sum(deg(v) choose 2) + {c5} * Sum(deg(v) choose 3)")
    print(f"The coefficients c_1, c_2, c_3, c_4, c_5 are:")
    print(f"{coefficients[0]}, {coefficients[1]}, {coefficients[2]}, {coefficients[3]}, {coefficients[4]}")

solve()