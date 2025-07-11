import math

def main():
    """
    This program calculates and prints the coefficients for the expression of the
    number of closed tree-like walks of length 6 in a simple graph X.

    The derivation is as follows:
    A closed tree-like walk of length 6 has an underlying graph structure that is a tree.
    Let the number of times edge i is traversed be n_i. Sum(n_i) = 6.
    For a closed walk on a tree, each n_i must be even. Let n_i = 2*k_i.
    So, Sum(k_i) = 3, with k_i >= 1.

    We analyze partitions of 3:
    1.  k_1 = 3: The tree has 1 edge (P_2), traversed 6 times.
        - For an edge u-v, walks are (u,v,u,v,u,v,u) and (v,u,v,u,v,u,v).
        - There are 2 such walks for each edge.
        - Number of edges = e. Total contribution: 2*e.
        - So, c_1 = 2.

    2.  k_1 = 2, k_2 = 1: The tree has 2 edges (P_3), traversed 4 and 2 times.
        - For a path u-v-w, let e1=(u,v) and e2=(v,w).
        - Walks with 4 traversals of e1 and 2 of e2: We found 6 such walks.
        - Walks with 2 traversals of e1 and 4 of e2: We found 6 such walks.
        - Total walks on a fixed P_3 is 12.
        - Number of P_3 subgraphs, centered at v, is C(deg(v), 2).
        - Total contribution: 12 * Sum(C(deg(v), 2)).
        - So, c_4 = 12.

    3.  k_1 = k_2 = k_3 = 1: The tree has 3 edges, each traversed twice.
        - Case a: Tree is a path P_4 (u-v-w-x).
          - We found 6 such walks on a fixed P_4.
          - Number of P_4 subgraphs = p. Total contribution: 6*p.
          - So, c_3 = 6.
        - Case b: Tree is a star K_1,3 (v central, connected to u,w,x).
          - We found 12 such walks on a fixed K_1,3.
          - Number of K_1,3 subgraphs, centered at v, is C(deg(v), 3).
          - Total contribution: 12 * Sum(C(deg(v), 3)).
          - So, c_5 = 12.

    4.  The k term (number of K_3 subgraphs).
        - A tree-like walk requires the graph of its edges to be a tree.
        - A K_3 (triangle) is a cycle, not a tree.
        - Therefore, there are no tree-like walks whose edge set forms a K_3.
        - So, c_2 = 0.
    """

    c1 = 2
    c2 = 0
    c3 = 6
    c4 = 12
    c5 = 12

    print("The coefficients for the expression are:")
    print(f"c_1 = {c1}")
    print(f"c_2 = {c2}")
    print(f"c_3 = {c3}")
    print(f"c_4 = {c4}")
    print(f"c_5 = {c5}")
    print("\nThe equation is:")
    print(f"Number of walks = {c1} * e + {c2} * k + {c3} * p + {c4} * Sum(C(deg(v), 2)) + {c5} * Sum(C(deg(v), 3))")

if __name__ == "__main__":
    main()