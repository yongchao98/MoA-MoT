def count_strata():
    """
    Calculates the number of codimension 2 boundary strata of the moduli space M_bar_{3,1}.

    This is equivalent to counting the number of topological types of stable curves of
    arithmetic genus 3 with 1 marked point and 2 nodes.

    The method is to enumerate all possible stable curve configurations.
    A configuration is defined by:
    1. The dual graph (number of components and how they are connected).
    2. The assignment of a geometric genus to each component.
    3. The placement of the marked point.

    The key constraints are:
    - The genus formula for a nodal curve: sum(g_i) = k_c (for g=3, delta=2)
    - The stability condition for each component: 2*g_i - 2 + n_i > 0
    """

    # List to hold the count for each distinct topological type. Each valid type adds 1.
    strata_counts = []
    
    # Case 1: k_c = 1 component (irreducible curve with two nodes).
    # The dual graph has 1 vertex and 2 loops. sum(g_i) = k_c => g_1 = 1.
    # The component has genus 1, 2 nodes (which count as 2 special points), and 1 marked point.
    # n_1 = 2 (nodes) + 1 (m.p.) = 3.
    # Stability: 2*1 - 2 + 3 = 3 > 0. This is a valid stratum.
    strata_counts.append(1)

    # Case 2: k_c = 2 components, connected by 2 nodes.
    # The dual graph has 2 vertices and 2 parallel edges. sum(g_i) = k_c => g_1 + g_2 = 2.
    # We analyze partitions of the genus sum.
    
    # Subcase 2a: Genus partition {2, 0}.
    # One component has g=2, the other g=0. The g=0 component needs at least 3 special points
    # to be stable. It has 2 node connections, so it must carry the marked point. This configuration is stable.
    strata_counts.append(1)

    # Subcase 2b: Genus partition {1, 1}.
    # Both components have g=1. The marked point on either component gives a symmetric configuration.
    # This configuration is stable.
    strata_counts.append(1)

    # Case 3: k_c = 3 components in a chain.
    # The dual graph is a line of 3 vertices and 2 edges. sum(g_i) = k_c => g_1 + g_2 + g_3 = 3.
    # Leaf components (ends of the chain) cannot be of genus 0 due to the stability condition.
    
    # Subcase 3a: Genus partition {2, 1, 0}.
    # The g=0 component must be in the middle of the chain and hold the marked point to be stable.
    # The two leaf components have genera 1 and 2. This gives one valid stratum type.
    strata_counts.append(1)
    
    # Subcase 3b: Genus partition {1, 1, 1}.
    # All components have g=1. All components are inherently stable. There are two distinct
    # topological positions for the single marked point.
    
    # 1. Marked point on a leaf component.
    strata_counts.append(1)
    
    # 2. Marked point on the central component.
    strata_counts.append(1)

    # The total number is the sum of counts for all valid topological types.
    total_strata = sum(strata_counts)
    
    # Output the result as an equation showing the contribution of each stratum type.
    equation_str = " + ".join(map(str, strata_counts))
    print(f"{equation_str} = {total_strata}")

count_strata()
<<<6>>>