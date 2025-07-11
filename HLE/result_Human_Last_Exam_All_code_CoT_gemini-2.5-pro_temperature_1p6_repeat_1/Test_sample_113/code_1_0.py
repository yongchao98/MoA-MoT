def is_stable(g, n_special_points):
    """Checks the stability condition for a curve component."""
    return 2 * g - 2 + n_special_points > 0

def solve():
    """
    Calculates the number of codimension 2 boundary strata of M_3,1.
    This is equivalent to counting stable dual graphs with 2 edges and 1 leg.
    """
    counts = []

    # Case 1: 1 vertex, 2 edges (must be loops)
    # g = h^1 + sum(g_v) => 3 = (2 - 1 + 1) + g_1 => g_1 = 1.
    # The single component has genus 1.
    # Number of special points on it: 1 leg + 2 * 2 connections from loops = 5.
    case1_count = 0
    if is_stable(g=1, n_special_points=5):
        case1_count = 1
    counts.append(case1_count)

    # Case 2: 2 vertices, 2 edges (must be parallel)
    # g = h^1 + sum(g_v) => 3 = (2 - 2 + 1) + (g_1 + g_2) => g_1 + g_2 = 2.
    # Genus partitions of 2: (2,0) and (1,1).
    case2_count = 0
    # Subcase 2a: Genera (2, 0).
    # Marked point on g=2 comp: v1(g=2,n=3), v2(g=0,n=2). v2 is unstable.
    # Marked point on g=0 comp: v1(g=2,n=2), v2(g=0,n=3). Both stable.
    if is_stable(g=2, n_special_points=2) and is_stable(g=0, n_special_points=3):
        case2_count += 1
    # Subcase 2b: Genera (1, 1).
    # Symmetric. Place marked point on v1. v1(g=1,n=3), v2(g=1,n=2). Both stable.
    if is_stable(g=1, n_special_points=3) and is_stable(g=1, n_special_points=2):
        case2_count += 1
    counts.append(case2_count)

    # Case 3: 3 vertices, 2 edges (must be a path v1-v2-v3)
    # g = h^1 + sum(g_v) => 3 = (2 - 3 + 1) + (g_1+g_2+g_3) => g_1+g_2+g_3 = 3.
    # Genus partitions of 3: (3,0,0), (2,1,0), (1,1,1).
    case3_count = 0
    # (3,0,0) yields no stable curves.
    # Subcase 3a: Genera (2,1,0). Only stable assignment is (1,0,2) or (2,0,1)
    # on the path, with marked point on the central g=0 component.
    # v_end(g=2,n=1), v_mid(g=0,n=3), v_end(g=1,n=1).
    if is_stable(g=2, n_special_points=1) and is_stable(g=0, n_special_points=3) and is_stable(g=1, n_special_points=1):
        case3_count += 1
    # Subcase 3b: Genera (1,1,1).
    # Point on an end vertex (symmetric): v_end(g=1,n=2), v_mid(g=1,n=2), v_end(g=1,n=1). All stable.
    if is_stable(g=1, n_special_points=2) and is_stable(g=1, n_special_points=2) and is_stable(g=1, n_special_points=1):
        case3_count += 1
    # Point on the middle vertex: v_end(g=1,n=1), v_mid(g=1,n=3), v_end(g=1,n=1). All stable.
    if is_stable(g=1, n_special_points=1) and is_stable(g=1, n_special_points=3) and is_stable(g=1, n_special_points=1):
        case3_count += 1
    counts.append(case3_count)

    total = sum(counts)
    equation = " + ".join(map(str, counts))
    print(f"The number of codimension 2 strata arise from three graph topologies.")
    print(f"Number of strata from a 1-vertex graph: {counts[0]}")
    print(f"Number of strata from a 2-vertex graph: {counts[1]}")
    print(f"Number of strata from a 3-vertex graph: {counts[2]}")
    print("\nThe total number of strata is the sum:")
    print(f"{equation} = {total}")

if __name__ == '__main__':
    solve()