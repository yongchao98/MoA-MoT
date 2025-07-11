import collections

def count_two_matchings(n, edges):
    """
    Counts the number of 2-matchings in a graph.

    A 2-matching is a set of 2 edges that do not share any vertices.
    This function implements an O(|V| + |E|) algorithm, which is subcubic.
    
    The formula used is:
    Count = (1/2) * sum_{(u,v) in E} (m - deg(u) - deg(v) + 1)
    where n = |V|, m = |E|, and deg(x) is the degree of vertex x.

    Args:
        n: The number of vertices in the graph (labeled 0 to n-1).
        edges: A list of tuples, where each tuple represents an edge.

    Returns:
        The number of 2-matchings in the graph.
    """
    if n == 0 or not edges:
        print("Graph is empty.")
        print("Number of 2-matchings: 0")
        return

    m = len(edges)
    degrees = collections.defaultdict(int)
    for u, v in edges:
        degrees[u] += 1
        degrees[v] += 1

    # This sum counts each pair of disjoint edges twice
    double_count_sum = 0
    for u, v in edges:
        # For an edge (u,v), the number of other edges not adjacent to it is
        # m - 1 (total edges minus itself)
        # - (degrees[u] - 1) (other edges connected to u)
        # - (degrees[v] - 1) (other edges connected to v)
        # This simplifies to m - degrees[u] - degrees[v] + 1
        disjoint_edges_count = m - degrees[u] - degrees[v] + 1
        double_count_sum += disjoint_edges_count
        
    # Each 2-matching {e1, e2} is counted once when we iterate on e1
    # and once when we iterate on e2. So we divide by 2.
    # Note: integer division // is used since the result must be an integer.
    count = double_count_sum // 2

    print(f"Graph properties:")
    print(f"Number of vertices n = {n}")
    print(f"Number of edges m = {m}")
    print("\nEquation for counting 2-matchings:")
    print("N_2 = (sum over each edge (u,v) of [m - deg(u) - deg(v) + 1]) / 2")
    print("\nCalculation steps:")
    print(f"The sum part of the equation evaluates to: {double_count_sum}")
    print(f"Final calculation: {double_count_sum} / 2 = {count}")
    print(f"\nResult: The number of 2-matchings is {count}.")


# Example usage: A cycle graph with 4 vertices (a square)
# V = {0, 1, 2, 3}, E = {(0,1), (1,2), (2,3), (3,0)}
# The 2-matchings are {(0,1), (2,3)} and {(1,2), (3,0)}. The count should be 2.
num_vertices = 4
edge_list = [(0, 1), (1, 2), (2, 3), (3, 0)]
count_two_matchings(num_vertices, edge_list)