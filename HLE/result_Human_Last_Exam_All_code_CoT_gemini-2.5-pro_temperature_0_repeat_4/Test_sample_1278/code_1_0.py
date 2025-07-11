import math

def main():
    """
    This script calculates the maximum number of Mexican standoffs (cycles)
    for the given problem conditions.

    The strategy is to construct a graph with 9 vertices and 16 edges that
    maximizes the number of cycles. The chosen structure consists of two dense
    subgraphs, G1 and G2, connected by two edges.

    - G1: A complete graph on 5 vertices with 1 edge removed (K_5-e).
          This has 5 vertices and 9 edges.
    - G2: A complete graph on 4 vertices with 1 edge removed (K_4-e).
          This has 4 vertices and 5 edges.
    - The remaining 2 edges connect the non-adjacent vertices in G1 to the
      non-adjacent vertices in G2.

    The total number of cycles is the sum of cycles within G1, cycles within G2,
    and the new cycles formed by combining paths from G1 and G2.
    """

    # --- Part 1: Cycles within G1 (K_5-e) ---

    # A full K_5 has:
    # 10 triangles (3-cycles)
    # 15 4-cycles
    # 12 5-cycles
    cycles_K5 = 10 + 15 + 12

    # Removing one edge from K_5 removes:
    # 3 triangles, 3 4-cycles, and 2 5-cycles.
    cycles_lost_in_K5_e = 3 + 3 + 2
    cycles_G1 = cycles_K5 - cycles_lost_in_K5_e

    # --- Part 2: Cycles within G2 (K_4-e) ---

    # A K_4-e graph has:
    # 2 triangles (3-cycles)
    # 1 4-cycle
    cycles_G2 = 2 + 1

    # --- Part 3: New cycles from connecting G1 and G2 ---

    # The number of new cycles is the product of the number of paths between
    # the connection points in each subgraph.

    # Number of paths between the two non-adjacent vertices in G1 (K_5-e):
    # Paths of length 2: 3
    # Paths of length 3: 6
    # Paths of length 4: 6
    paths_G1 = 3 + 6 + 6

    # Number of paths between the two non-adjacent vertices in G2 (K_4-e):
    # Paths of length 2: 2
    # Paths of length 3: 2
    paths_G2 = 2 + 2

    new_cycles = paths_G1 * paths_G2

    # --- Final Calculation ---
    total_cycles = cycles_G1 + cycles_G2 + new_cycles

    print("The maximum number of Mexican standoffs is found by summing the cycles from three sources:")
    print("1. Cycles within the first subgraph (a K_5-e): {}".format(cycles_G1))
    print("2. Cycles within the second subgraph (a K_4-e): {}".format(cycles_G2))
    print("3. New cycles formed by connecting the subgraphs: {} paths in G1 * {} paths in G2 = {}".format(paths_G1, paths_G2, new_cycles))
    print("\nFinal equation:")
    print("{} + {} + {} * {} = {}".format(cycles_G1, cycles_G2, paths_G1, paths_G2, total_cycles))

if __name__ == "__main__":
    main()