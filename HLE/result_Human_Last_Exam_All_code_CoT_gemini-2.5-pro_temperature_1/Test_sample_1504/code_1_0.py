import sys

def solve_graph_theory_problem():
    """
    Solves the given Turan-type extremal graph theory problem.

    The reasoning is as follows:

    (a) True. If G is not a union of K_2's, it contains a connected component with at least 3 vertices.
        - Lower bound: A matching on n vertices has floor(n/2) edges, is G-free, and is K_{1,t}-induced-free for t>=2. This gives an Omega(n) lower bound.
        - Upper bound: Let H be a (G, K_{1,t}-ind)-free graph. For any vertex v, its neighborhood N(v) must not contain an independent set of size t. By Ramsey's Theorem, if the degree of v is at least R(|V(G)|, t), then N(v) must contain a clique of size |V(G)|. This clique, along with v, forms a subgraph that contains G. Thus, the maximum degree of H is bounded by the constant R(|V(G)|, t) - 1. A graph with constantly bounded maximum degree has O(n) edges.
        - Since the function is both Omega(n) and O(n), it is Theta(n).

    (b) True. If a graph H is sK_2-free, its maximum matching size is at most s-1.
        - Let H' be the subgraph of H containing all non-isolated vertices. Let M be a maximum matching in H'; |V(M)| <= 2(s-1). Let I = V(H') - V(M); I is an independent set.
        - Any vertex x in V(M) can be adjacent to fewer than t vertices in I, otherwise, those neighbors would form an independent set of size t in the neighborhood of x, which is forbidden.
        - This implies the total number of vertices in I is also bounded by a function of s and t.
        - Therefore, the total number of non-isolated vertices in H is bounded by a constant depending on s and t. The number of edges is thus also bounded by a constant. For n larger than this constant, the maximum number of edges does not grow, making it independent of n.

    (c) The tight upper bound is known to be the maximum of the number of edges in two specific constructions: K_{2s-1} and K_{s-1} joined with an independent set of size t-1 (K_{s-1} v I_{t-1}). The latter's edge count provides a general expression involving both s and t.
        - The number of edges in K_{s-1} v I_{t-1} is the sum of edges in K_{s-1} and the edges between the two sets.
        - Edges = C(s-1, 2) + (s-1)(t-1)
        - Simplifying this expression:
          = (s-1)(s-2)/2 + (s-1)(t-1)
          = (s-1) * [(s-2)/2 + (t-1)]
          = (s-1) * [(s-2 + 2t - 2)/2]
          = (s-1)(s+2t-4)/2
    """
    
    # Part (a)
    answer_a = "True"
    
    # Part (b)
    answer_b = "True"
    
    # Part (c)
    # The expression is derived from the number of edges in the graph K_{s-1} v I_{t-1}.
    # |E| = C(s-1, 2) + (s-1)(t-1) = (s-1)(s-2)/2 + (s-1)(t-1)
    # This simplifies to (s-1)(s+2t-4)/2
    expression_c = "(s-1)(s+2t-4)/2"
    
    # Printing the final answer in the required format
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"
    print(final_answer)

solve_graph_theory_problem()