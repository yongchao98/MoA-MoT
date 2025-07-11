def solve():
    """
    This function determines the coefficients for the number of closed tree-like walks of length 6.

    The number of closed tree-like walks of length 6 in a simple graph X is given by:
    c_1 * e + c_2 * k + c_3 * p + c_4 * sum(deg(v) choose 2) + c_5 * sum(deg(v) choose 3)

    where:
    e: number of edges
    k: number of K_3 subgraphs (triangles)
    p: number of P_4 subgraphs (paths of length 3)
    deg(v): degree of vertex v

    Our analysis proceeds by classifying tree-like walks based on the structure of the tree formed by their edges.
    A closed walk of length 6 on a tree must traverse each edge an even number of times.

    1.  Tree is K_2 (1 edge): Edge traversed 6 times. 2 walks per edge. Contribution: 2*e. So, c_1 = 2.
    2.  Tree is P_3 (2 edges): One edge traversed 4 times, the other 2 times. 12 walks per P_3. Contribution: 12 * sum(deg(v) choose 2). So, c_4 = 12.
    3.  Tree is P_4 (3 edges): Each edge traversed twice. 6 walks per P_4. Contribution: 6*p. So, c_3 = 6.
    4.  Tree is K_1,3 (3 edges): Each edge traversed twice. 12 walks per K_1,3. Contribution: 12 * sum(deg(v) choose 3). So, c_5 = 12.
    5.  K_3 term: A tree-like walk cannot have an underlying edge set that contains a cycle. Therefore, walks whose fundamental structure is a K_3 are not tree-like. Hence, c_2 = 0.
    """
    
    c1 = 2
    c2 = 0
    c3 = 6
    c4 = 12
    c5 = 12
    
    print(c1)
    print(c2)
    print(c3)
    print(c4)
    print(c5)

solve()