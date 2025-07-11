import collections
import itertools

memo = {}

def get_g_recursive(graph_edges):
    """
    Computes the Grundy value (g-value) of a game state represented by a graph
    using the recursive mex definition.
    A graph is represented by a frozenset of its edges.
    """
    if not graph_edges:
        return 0
    
    # Use a canonical (sorted tuple) representation for memoization keys
    graph_key = tuple(sorted(list(edge) for edge in graph_edges))
    if graph_key in memo:
        return memo[graph_key]

    adj = collections.defaultdict(set)
    vertices = set()
    for u, v in graph_edges:
        adj[u].add(v)
        adj[v].add(u)
        vertices.add(u)
        vertices.add(v)

    reachable_g_values = set()
    
    # A move is to pick a vertex and remove a non-empty subset of its incident edges
    for v in vertices:
        incident_edges = [frozenset([v, neighbor]) for neighbor in adj[v]]
        
        # Iterate through all non-empty subsets of incident edges
        for i in range(1, 1 << len(incident_edges)):
            edges_to_remove = frozenset(
                itertools.compress(incident_edges, ((i >> j) & 1 for j in range(len(incident_edges))))
            )
            
            next_graph_edges = graph_edges - edges_to_remove
            
            # The g-value of the resulting state is the nim-sum of its components' g-values
            g_val = get_g_optimized(next_graph_edges)
            reachable_g_values.add(g_val)

    # Compute mex (Minimum Excluded value)
    mex = 0
    while mex in reachable_g_values:
        mex += 1
    
    memo[graph_key] = mex
    return mex

def get_g_optimized(graph_edges):
    """
    Computes the g-value of a graph by decomposing it into connected components
    and using known results for simpler structures like trees.
    """
    if not graph_edges:
        return 0

    adj = collections.defaultdict(set)
    vertices = set()
    for u, v in graph_edges:
        adj[u].add(v)
        adj[v].add(u)
        vertices.add(u)
        vertices.add(v)

    visited = set()
    nim_sum = 0
    
    # Find connected components and compute nim-sum of their g-values
    for start_node in vertices:
        if start_node not in visited:
            component_nodes = set()
            component_edges = set()
            q = collections.deque([start_node])
            visited.add(start_node)
            
            # BFS to find all nodes and edges in the component
            while q:
                curr = q.popleft()
                component_nodes.add(curr)
                for neighbor in adj[curr]:
                    component_edges.add(frozenset([curr, neighbor]))
                    if neighbor not in visited:
                        visited.add(neighbor)
                        q.append(neighbor)
            
            # Optimization: for trees, g-value is the number of edges
            if len(component_edges) == len(component_nodes) - 1:
                nim_sum ^= len(component_edges)
            else:
                nim_sum ^= get_g_recursive(frozenset(component_edges))

    return nim_sum
    
def solve(n, m):
    """
    Calculates f(n, m) by enumerating all possible matrices, computing the
    g-value for each, and checking if the winning probability > 0.5.
    """
    winning_configs = 0
    total_configs = 1 << (n * m)

    # Iterate through all 2^(n*m) binary matrices
    for i in range(total_configs):
        graph_edges = set()
        # Build the bipartite graph for the current matrix
        for r in range(n):
            for c in range(m):
                if (i >> (r * m + c)) & 1:
                    # Row vertices are 'r0', 'r1', ...; Column vertices are 'c0', 'c1', ...
                    graph_edges.add(frozenset([f'r{r}', f'c{c}']))
        
        # Calculate the g-value for the current game state
        g_value = get_g_optimized(frozenset(graph_edges))

        if g_value > 0:
            winning_configs += 1

    probability = winning_configs / total_configs if total_configs > 0 else 0
    
    print(f"Analysis for a {n}x{m} matrix:")
    print(f"Total configurations: {total_configs}")
    print(f"Winning configurations (g > 0): {winning_configs}")
    print(f"Losing configurations (g = 0): {total_configs - winning_configs}")
    print(f"Probability of 1st player winning: {probability:.4f}")

    if probability > 0.5:
        print("Probability is > 0.5, so the first player has the advantage.")
        return 1
    else:
        print("Probability is not > 0.5, so the first player does not have the advantage.")
        return 0

# Example usage for a small case.
# Note: This will be very slow for n*m > 8
# solve(2, 2) 

# The problem asks for the complexity, not the result of a specific run.
# Based on the analysis, the complexity is determined by the enumeration
# of all matrices and the exponential-time calculation of the g-value.
print("The computational complexity of the function f(n, m) is determined by an algorithm that enumerates all 2^(n*m) possible game states and computes the Grundy value for each.")
print("The Grundy value calculation is itself exponential in the worst case.")
print("This results in a total complexity that is exponential in the size of the matrix.")
