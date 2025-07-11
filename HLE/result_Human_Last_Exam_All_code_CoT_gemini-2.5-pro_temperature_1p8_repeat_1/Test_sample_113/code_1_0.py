def count_codim_2_strata_M31():
    """
    Calculates the number of codimension 2 boundary strata of M_bar(3,1)
    by enumerating all possible decorated dual graphs.
    """
    print("Enumerating codimension 2 boundary strata of M_bar(3,1):")
    print("This corresponds to stable curves with 2 nodes.")
    print("-" * 50)

    total_strata = 0
    strata_counts = []

    # Case 1: Irreducible curve (1 component)
    # The dual graph has |V|=1 vertex and |E|=2 edges (loops).
    # Genus constraint: g_1 = |V| = 1.
    # The component must have genus 1.
    # Stability: g=1, 1 marked point, 2 self-nodes.
    # Number of special points n_1 = 1 (marked point) + 2*2 (for two loops) = 5.
    # Stability check: 2*1 - 2 + 5 = 5 > 0. This is a valid configuration.
    num_strata_v1 = 1
    print(f"Case 1: Irreducible curves (1 vertex, 2 loops)")
    print(f"  - Genus assignment: (g=1). Unique and stable.")
    print(f"  - Number of strata: {num_strata_v1}")
    total_strata += num_strata_v1
    strata_counts.append(num_strata_v1)
    print("-" * 50)

    # Case 2: Curve with 2 components
    # The dual graph has |V|=2 vertices and |E|=2 edges.
    # Genus constraint: g_1 + g_2 = |V| = 2.
    # Two possible graph structures: two parallel edges, or one edge and a loop.

    # Subcase 2a: Two parallel edges connecting the two vertices.
    # One vertex has the marked point. Let this be v1. The other, v2, does not.
    # n1 = 2 (nodes) + 1 (mark) = 3. Stability: 2g1-2+3 > 0 -> always stable.
    # n2 = 2 (nodes). Stability: 2g2-2+2 > 0 -> g2 > 0.
    # We need to find pairs (g1, g2) with g1+g2=2, g1>=0, g2>=1.
    # (g1,g2) = (1,1): Valid. One stratum.
    # (g1,g2) = (0,2): Valid. One stratum.
    num_strata_v2a = 2
    print(f"Case 2: Reducible curves (2 vertices, g1+g2=2)")
    print(f"  Subcase 2a: Two parallel edges (v1 -- v2)")
    print(f"    - (g1=1, g2=1) with marked point on one: 1 stratum")
    print(f"    - (g1=0, g2=2) with marked point on g=0 component (for stability of other): 1 stratum")
    print(f"  - Number of strata for this graph: {num_strata_v2a}")
    total_strata += num_strata_v2a
    strata_counts.append(num_strata_v2a)
    
    # Subcase 2b: One edge between vertices, one loop on one vertex (say v1).
    # The vertices are distinct: v1 has a loop, v2 does not.
    num_strata_v2b = 0
    #   i) Marked point on v1 (loopy vertex)
    #      n1 = 2(loop) + 1(edge) + 1(mark) = 4. Stability: 2g1-2+4>0 -> g1>=0.
    #      n2 = 1(edge). Stability: 2g2-2+1>0 -> g2>=1.
    #      (g1,g2)=(1,1) -> 1 stratum. (g1,g2)=(0,2) -> 1 stratum.
    num_strata_v2b += 2
    #   ii) Marked point on v2 (non-loopy vertex)
    #       n1 = 2(loop) + 1(edge) = 3. Stability: 2g1-2+3>0 -> g1>=0.
    #       n2 = 1(edge) + 1(mark) = 2. Stability: 2g2-2+2>0 -> g2>=1.
    #       (g1,g2)=(1,1) -> 1 stratum. (g1,g2)=(0,2) -> 1 stratum.
    num_strata_v2b += 2
    print(f"  Subcase 2b: One edge and one loop on v1 (v1 O-- v2)")
    print(f"    - Marked point on v1 (loopy): (g1,g2)=(1,1) or (0,2). 2 strata.")
    print(f"    - Marked point on v2 (non-loopy): (g1,g2)=(1,1) or (0,2). 2 strata.")
    print(f"  - Number of strata for this graph: {num_strata_v2b}")
    total_strata += num_strata_v2b
    strata_counts.append(num_strata_v2b)
    print("-" * 50)
    
    # Case 3: Curve with 3 components
    # The dual graph has |V|=3 vertices and |E|=2 edges. It must be a chain.
    # Genus constraint: g_1 + g_2 + g_3 = |V| = 3.
    # Let the chain be v1 -- v2 -- v3.
    num_strata_v3 = 0
    #   i) Marked point on an end vertex (e.g., v1).
    #      n1=2, n2=2, n3=1.
    #      Stability requires: g1>=1, g2>=1, g3>=1.
    #      The only partition is (g1,g2,g3) = (1,1,1). The ends are symmetric.
    num_strata_v3 += 1
    #   ii) Marked point on the middle vertex (v2).
    #       n1=1, n2=3, n3=1.
    #       Stability requires: g1>=1, g3>=1, g2>=0.
    #       Possible partitions (g1,g2,g3):
    #       - (1,1,1): Symmetric ends. 1 stratum.
    #       - (1,0,2): Asymmetric ends. 1 stratum.
    #       - (2,0,1): Distinct from previous as the ends are ordered. 1 stratum.
    num_strata_v3 += 3
    print(f"Case 3: Reducible curves (3 vertices, g1+g2+g3=3)")
    print(f"  - The graph is a chain v1 -- v2 -- v3.")
    print(f"    - Marked point on end (v1 or v3): (1,1,1). 1 stratum.")
    print(f"    - Marked point on middle (v2): (1,1,1), (1,0,2), (2,0,1). 3 strata.")
    print(f"  - Number of strata for this graph: {num_strata_v3}")
    total_strata += num_strata_v3
    strata_counts.append(num_strata_v3)
    print("-" * 50)
    
    print("Total number of strata is the sum of all cases.")
    equation_parts = [str(c) for c in strata_counts]
    print(f"{' + '.join(equation_parts)} = {total_strata}")

    return total_strata

if __name__ == '__main__':
    final_answer = count_codim_2_strata_M31()
    # The final answer needs to be enclosed in <<<>>>
    # The problem requires the code block to be self-contained, and the final output format.
    # As an AI model, I am adding the <<<...>>> format after the code block.
    # To run this code and see the logic, you can execute it as a python script.
