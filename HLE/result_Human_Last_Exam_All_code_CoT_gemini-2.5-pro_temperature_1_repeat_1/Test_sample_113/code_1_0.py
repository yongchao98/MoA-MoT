import math

def count_codim_2_strata_M_3_1():
    """
    Calculates the number of codimension 2 boundary strata of the moduli space M_bar_{3,1}.

    This corresponds to counting the number of stable dual graphs with g=3, n=1, and 2 edges.
    The method is to analyze the three possible connected graph structures with 2 edges.
    """
    print("This program calculates the number of codimension 2 boundary strata of M_bar_{3,1}.")
    print("This is equivalent to counting the number of topological types of stable curves of genus 3 with 1 marked point and 2 nodes.")
    print("-" * 70)

    # --- Case 1: Irreducible curve with two nodes ---
    # Dual graph: 1 vertex, 2 self-loops.
    # h^1(Gamma) = 2.
    # Genus sum: g_v = g - h^1 = 3 - 2 = 1.
    # The single component must have genus 1.
    # Special points on the component: 1 marked point + 2 nodes = 3.
    # Stability: 2*g_v - 2 + num_special_points > 0
    #           2*1 - 2 + 3 = 3 > 0. This is stable.
    count1 = 1
    print("Case 1: Irreducible curve with two nodes.")
    print("This corresponds to a single component of geometric genus 1, with two nodes and one marked point.")
    print(f"Number of strata in this case: {count1}")
    print("-" * 70)


    # --- Case 2: Two components joined at two points ---
    # Dual graph: 2 vertices, 2 edges between them.
    # h^1(Gamma) = 1.
    # Genus sum: g1 + g2 = g - h^1 = 3 - 1 = 2.
    # Each component is connected to 2 nodes.
    
    # We enumerate distinct strata based on component genera (g1, g2) and marked point location.
    # A stratum is defined by the set of its components { (g, num_points, num_nodes) }.
    strata_case2 = set()

    # Iterate over partitions of 2 into g1+g2
    for g1 in range(3):
        g2 = 2 - g1
        
        # Possibility A: Marked point is on component 1
        # Comp 1: genus g1, 1 marked point, 2 nodes. Stability: 2*g1 - 2 + (1+2) > 0 -> 2*g1 + 1 > 0 (always true)
        # Comp 2: genus g2, 0 marked points, 2 nodes. Stability: 2*g2 - 2 + 2 > 0 -> g2 > 0
        if g2 > 0:
            # Create a canonical representation of the stratum: frozenset of component tuples (genus, num_points)
            stratum = frozenset([(g1, 1), (g2, 0)])
            strata_case2.add(stratum)

        # Possibility B: Marked point is on component 2
        # Similar check, but now g1 > 0 is required for stability of the component without the marked point.
        if g1 > 0:
            stratum = frozenset([(g1, 0), (g2, 1)])
            strata_case2.add(stratum)
            
    count2 = len(strata_case2)
    print("Case 2: A curve with two components joined at two nodes.")
    print("The sum of the genera of the components must be 2. We check all stable configurations:")
    print(" - A genus 1 curve with the point, attached to a genus 1 curve.")
    print(" - A rational (g=0) curve with the point, attached to a genus 2 curve.")
    print(" - A genus 2 curve with the point, attached to a rational (g=0) curve.")
    print(f"Number of strata in this case: {count2}")
    print("-" * 70)


    # --- Case 3: A chain of three components ---
    # Dual graph: 3 vertices, 2 edges in a line.
    # h^1(Gamma) = 0.
    # Genus sum: g1 + g2 + g3 = g - h^1 = 3 - 0 = 3.
    # The end components (v1, v3) have 1 node. The middle component (v2) has 2 nodes.
    
    strata_case3 = set()

    # Iterate over partitions of 3 into g1+g2+g3
    for g1 in range(4):
        for g2 in range(4 - g1):
            g3 = 3 - g1 - g2
            
            # Possibility A: Marked point on an end component (v1 by symmetry)
            # v1: g1, 1 point, 1 node. Stability: 2*g1 - 2 + (1+1) > 0 -> g1 > 0
            # v2: g2, 0 points, 2 nodes. Stability: 2*g2 - 2 + 2 > 0 -> g2 > 0
            # v3: g3, 0 points, 1 node. Stability: 2*g3 - 2 + 1 > 0 -> g3 > 0
            if g1 > 0 and g2 > 0 and g3 > 0:
                # Canonical representation: ('end_point', frozenset of end genera, middle genus)
                # Since the end components are symmetric, we use a frozenset for their genera.
                stratum = ('end_point', frozenset([g1, g3]), g2)
                strata_case3.add(stratum)

            # Possibility B: Marked point on the middle component (v2)
            # v1: g1, 0 points, 1 node. Stability: 2*g1 - 2 + 1 > 0 -> g1 > 0
            # v2: g2, 1 point, 2 nodes. Stability: 2*g2 - 2 + (1+2) > 0 -> 2*g2 + 1 > 0 (always true)
            # v3: g3, 0 points, 1 node. Stability: 2*g3 - 2 + 1 > 0 -> g3 > 0
            if g1 > 0 and g3 > 0:
                # Canonical representation: ('mid_point', sorted tuple of end genera, middle genus)
                # Here g1 and g3 are distinguishable if different, but the overall curve is symmetric
                # if we swap them. So we sort them to get a canonical key.
                stratum = ('mid_point', tuple(sorted((g1, g3))), g2)
                strata_case3.add(stratum)

    count3 = len(strata_case3)
    print("Case 3: A curve made of three components connected in a chain.")
    print("The sum of the genera of the components must be 3. We check all stable configurations:")
    print(" - Three genus 1 curves in a chain, point on an end component.")
    print(" - Three genus 1 curves in a chain, point on the central component.")
    print(" - A genus 1, a rational (g=0), and a genus 2 curve in a chain, point on the rational component.")
    print(f"Number of strata in this case: {count3}")
    print("-" * 70)

    # --- Final Summation ---
    total = count1 + count2 + count3
    print("Total number of codimension 2 boundary strata is the sum of the counts from all cases.")
    print(f"Final Calculation: {count1} + {count2} + {count3} = {total}")
    
    return total

if __name__ == '__main__':
    final_answer = count_codim_2_strata_M_3_1()
    # The final answer is printed within the function.
    # To conform to the output format, we also print it here.
    # print(f"<<<{final_answer}>>>")