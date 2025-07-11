import networkx as nx
from itertools import combinations

def has_dominating_set(G, l):
    """
    Checks if a graph G has a dominating set of size l.
    This is a naive, brute-force implementation.
    """
    nodes = list(G.nodes())
    n = len(nodes)
    if l <= 0:
        return 0

    for combo in combinations(nodes, l):
        dominated_nodes = set(combo)
        for u in combo:
            dominated_nodes.update(G.neighbors(u))
        
        if len(dominated_nodes) == n:
            return 1 # Found a dominating set
    return 0 # No dominating set of size l found

def count_independent_sets(G, l):
    """
    Counts the number of independent sets of size l in a graph G.
    This is a naive, brute-force implementation.
    """
    nodes = list(G.nodes())
    count = 0
    if l <= 0:
        return 0
        
    for combo in combinations(nodes, l):
        is_independent = True
        # Check all pairs of vertices in the combination for an edge
        for u, v in combinations(combo, 2):
            if G.has_edge(u, v):
                is_independent = False
                break
        if is_independent:
            count += 1
    return count

# Create a simple graph to demonstrate the functions
# A cycle graph on 6 vertices: 0-1-2-3-4-5-0
G = nx.cycle_graph(6)
l_dom = 2 # A dominating set of size 2, e.g., {0, 3}
l_ind = 3 # An independent set of size 3, e.g., {0, 2, 4}

dom_set_result = has_dominating_set(G, l_dom)
ind_set_count = count_independent_sets(G, l_ind)

print(f"Problem 1: DomSet")
print(f"Graph has {G.number_of_nodes()} vertices.")
print(f"Does the graph have a dominating set of size {l_dom}?")
print(f"Output: {dom_set_result}")

print(f"\nProblem 2: #IndSet")
print(f"Graph has {G.number_of_nodes()} vertices.")
print(f"What is the number of independent sets of size {l_ind}?")
print(f"Output: {ind_set_count}")
