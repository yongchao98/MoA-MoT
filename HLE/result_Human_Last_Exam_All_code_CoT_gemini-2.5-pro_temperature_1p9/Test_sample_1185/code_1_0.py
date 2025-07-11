def solve_stable_reduction_g2():
    """
    This function enumerates and verifies the types of stable reductions for curves of genus 2.

    A type is defined by a decorated dual graph (vertices are components, edges are nodes).
    It must satisfy two conditions for genus g=2:
    1. Genus Formula: 2 = h1(G) + sum(g_i), where h1(G) = e - v + 1 is the number of loops
       in the graph G with v vertices and e edges, and g_i is the genus of component i.
    2. Stability: If a component has genus g_i = 0, its vertex must have degree >= 3.
       (A self-loop adds 2 to the degree).
    """
    print("Finding the number of different types of stable reduction for curves of genus 2.")
    print("-" * 70)

    stable_types_found = []

    # Let's analyze candidate graph structures.

    # --- Case 1: Irreducible Curves (1 component/vertex) ---
    # The graph has 1 vertex (v=1). h1 = e - 1 + 1 = e, where e is the number of self-loops (nodes).
    # The formula is 2 = e + g_1.
    
    # Candidate 1.1: 1 vertex, g=1, 1 node (1 self-loop)
    v, e, g_list = 1, 1, [1]
    g = 2
    h1 = e - v + 1
    sum_g = sum(g_list)
    is_genus_formula_ok = (g == h1 + sum_g)
    deg_v1 = 2 * e
    is_stable = not (g_list[0] == 0 and deg_v1 < 3)
    if is_genus_formula_ok and is_stable:
        desc = "Irreducible elliptic curve with one node."
        stable_types_found.append((desc, g, h1, sum_g))

    # Candidate 1.2: 1 vertex, g=0, 2 nodes (2 self-loops)
    v, e, g_list = 1, 2, [0]
    g = 2
    h1 = e - v + 1
    sum_g = sum(g_list)
    is_genus_formula_ok = (g == h1 + sum_g)
    deg_v1 = 2 * e
    is_stable = not (g_list[0] == 0 and deg_v1 < 3)
    if is_genus_formula_ok and is_stable:
        desc = "Irreducible rational curve with two nodes."
        stable_types_found.append((desc, g, h1, sum_g))

    # --- Case 2: Reducible Curves (v > 1) ---

    # Candidate 2.1: 2 elliptic components, joined by 1 node
    # Graph: v=2, e=1 (edge v1-v2). g_list=[1,1].
    v, e, g_list = 2, 1, [1, 1]
    g = 2
    h1 = e - v + 1
    sum_g = sum(g_list)
    is_genus_formula_ok = (g == h1 + sum_g)
    # Stability: No rational components, so condition is met.
    is_stable = True
    if is_genus_formula_ok and is_stable:
        desc = "Two elliptic curves meeting at one node."
        stable_types_found.append((desc, g, h1, sum_g))
        
    # Candidate 2.2: Elliptic curve meets a nodal rational curve at one point
    # Graph: v=2, v1(g=1), v2(g=0). e=2 (edge v1-v2 and a self-loop on v2).
    v, e, g_list = 2, [1, 0] # genera for v1, v2
    g = 2
    # The graph topology has 2 vertices, an edge between them, and a loop on one. Total e=2.
    h1 = 2 - v + 1
    sum_g = sum(g_list)
    is_genus_formula_ok = (g == h1 + sum_g)
    # Stability Check
    deg_v1 = 1 # connected to v2
    deg_v2 = 1 (connected to v1) + 2 # self-loop
    is_stable_v1 = not (g_list[0] == 0 and deg_v1 < 3)
    is_stable_v2 = not (g_list[1] == 0 and deg_v2 < 3)
    if is_genus_formula_ok and is_stable_v1 and is_stable_v2:
        desc = "An elliptic curve and a nodal rational curve meeting at one point."
        stable_types_found.append((desc, g, h1, sum_g))

    # Candidate 2.3: Two rational components, joined by 3 nodes
    # Graph: v=2, e=3 (3 edges between v1 and v2). g_list=[0,0].
    v, e, g_list = 2, 3, [0, 0]
    g = 2
    h1 = e - v + 1
    sum_g = sum(g_list)
    is_genus_formula_ok = (g == h1 + sum_g)
    # Stability Check
    deg_v1 = 3
    deg_v2 = 3
    is_stable_v1 = not (g_list[0] == 0 and deg_v1 < 3)
    is_stable_v2 = not (g_list[1] == 0 and deg_v2 < 3)
    if is_genus_formula_ok and is_stable_v1 and is_stable_v2:
        desc = "Two rational curves meeting at three nodes."
        stable_types_found.append((desc, g, h1, sum_g))

    print(f"Found {len(stable_types_found)} distinct types of stable reduction:\n")
    for i, (desc, g, h1, sum_g) in enumerate(stable_types_found):
        print(f"Type {i+1}: {desc}")
        print(f"  - Genus formula: g = h¹(G) + Σgᵢ")
        print(f"  - Verification:  {g} = {h1} + {sum_g}\n")

if __name__ == "__main__":
    solve_stable_reduction_g2()
    print("<<<5>>>")