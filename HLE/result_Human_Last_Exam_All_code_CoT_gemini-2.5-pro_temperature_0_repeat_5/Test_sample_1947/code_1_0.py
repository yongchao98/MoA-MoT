import math

def get_coefficients():
    """
    Calculates the coefficients c_1 through c_5 based on combinatorial counting of closed tree-like walks of length 6.
    """

    # c_1: Walks on a single edge (K_2)
    # The underlying tree is a single edge (u,v). The walk length is 6.
    # This means the edge must be traversed 3 times in each direction.
    # If starting at u, the only possible walk is u-v-u-v-u-v-u. (1 walk)
    # If starting at v, the only possible walk is v-u-v-u-v-u-v. (1 walk)
    # Total walks per edge = 2.
    # The term is c_1 * e, where e is the number of edges. So, c_1 = 2.
    c1 = 2

    # c_2: Walks related to triangles (K_3)
    # A "tree-like" walk is one where the set of edges traversed forms a tree.
    # A K_3 (triangle) is not a tree as it contains a cycle.
    # Therefore, by this definition, there are no tree-like walks whose underlying graph is a K_3.
    # This implies the coefficient c_2 is 0.
    c2 = 0

    # c_3: Walks on a path of 3 edges (P_4)
    # The tree is a path a-b-c-d. Each of the 3 edges is traversed once in each direction.
    # Start at 'a': The only walk is a->b->c->d->c->b->a. (1 walk)
    # Start at 'd': By symmetry, 1 walk.
    # Start at 'b': The walk consists of a round trip to 'a' (b->a->b) and a round trip over the rest of the path (b->c->d->c->b).
    # These two parts can be done in any order, giving 2 walks:
    # 1. b->a->b -> c->d->c->b
    # 2. b->c->d->c->b -> a->b
    # Start at 'c': By symmetry with 'b', 2 walks.
    # Total walks per P_4 = 1 + 1 + 2 + 2 = 6.
    # The term is c_3 * p, where p is the number of P_4 subgraphs. So, c_3 = 6.
    c3 = 6

    # c_4: Walks on a path of 2 edges (P_3)
    # The tree is a path u-v-w. It has 2 edges, e1=(u,v) and e2=(v,w).
    # The sum of half-traversals is 3, so one edge is traversed 2*1=2 times and the other 2*2=4 times.
    # Case A: e1 traversed 4 times, e2 traversed 2 times.
    #   Start at u: Walk is u->v...v->u. Inner walk at v of length 4 uses one v<->u trip and one v<->w trip.
    #     Two orders for inner trips: v->u->v->w->v and v->w->v->u->v. This gives 2 full walks.
    #   Start at w: Walk is w->v...v->w. Inner walk at v of length 4 must be two v<->u trips: v->u->v->u->v. (1 walk)
    #   Start at v: Concatenation of two v<->u trips and one v<->w trip. Permutations of (Tu, Tu, Tw) is 3!/2! = 3 walks.
    #   Total for Case A = 2 + 1 + 3 = 6 walks.
    # Case B: e1 traversed 2 times, e2 traversed 4 times. By symmetry, 6 walks.
    # Total walks per P_3 = 6 + 6 = 12.
    # The term is c_4 * sum(deg(v) choose 2). So, c_4 = 12.
    c4 = 12

    # c_5: Walks on a star graph with 3 edges (K_1,3)
    # The tree has a central vertex v0 and 3 leaves v1, v2, v3. Each of the 3 edges is traversed once in each direction.
    # Start at center v0: The walk consists of three round trips v0<->vi. The order can be permuted in 3! = 6 ways. (6 walks)
    # Start at a leaf v1: Walk is v1->v0...v0->v1. The inner walk at v0 of length 4 consists of round trips to v2 and v3.
    #   The two trips can be ordered in 2! = 2 ways.
    #   There are 3 leaves, so 3 * 2 = 6 walks starting from leaves.
    # Total walks per K_1,3 = 6 + 6 = 12.
    # The term is c_5 * sum(deg(v) choose 3). So, c_5 = 12.
    c5 = 12
    
    return c1, c2, c3, c4, c5

def main():
    """
    Main function to get and print the coefficients.
    """
    c1, c2, c3, c4, c5 = get_coefficients()
    
    print(f"The coefficient c_1 is: {c1}")
    print(f"The coefficient c_2 is: {c2}")
    print(f"The coefficient c_3 is: {c3}")
    print(f"The coefficient c_4 is: {c4}")
    print(f"The coefficient c_5 is: {c5}")
    
    print("\nThe expression is:")
    print(f"{c1} * e + {c2} * k + {c3} * p + {c4} * sum(deg(v) choose 2) + {c5} * sum(deg(v) choose 3)")

if __name__ == "__main__":
    main()
