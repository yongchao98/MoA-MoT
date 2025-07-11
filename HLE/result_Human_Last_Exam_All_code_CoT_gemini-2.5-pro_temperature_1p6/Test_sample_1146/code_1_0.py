import math
from collections import deque

def get_units(n):
    """Computes the set of units in Z_n."""
    return {i for i in range(1, n) if math.gcd(i, n) == 1}

def build_associate_graph(n):
    """Builds the associate graph AG(Z_n)."""
    if n <= 1:
        return {}
    
    vertices = list(range(1, n))
    adj = {v: [] for v in vertices}
    units = get_units(n)

    # Check all pairs of distinct vertices for the associate relationship.
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            a = vertices[i]
            b = vertices[j]
            
            # Check if a and b are associates: a = bu (mod n) for some unit u.
            # This is equivalent to checking if a is in the orbit of b under
            # the multiplicative action of the group of units.
            is_associate = False
            for u in units:
                if (b * u) % n == a:
                    is_associate = True
                    break
            
            if is_associate:
                adj[a].append(b)
                adj[b].append(a)
                
    return adj

def is_cycle_graph(adj):
    """Checks if a graph is a single cycle."""
    n_vertices = len(adj)
    
    # A cycle must have at least 3 vertices.
    if n_vertices < 3:
        return False

    # Check 1: Every vertex in a cycle must have a degree of 2.
    for v in adj:
        if len(adj[v]) != 2:
            return False

    # Check 2: The graph must be connected.
    # A graph where every vertex has degree 2 is a disjoint union of cycles.
    # If it is also connected, it must be a single cycle.
    q = deque([next(iter(adj))])  # Start BFS from an arbitrary vertex.
    visited = {q[0]}
    
    while q:
        u = q.popleft()
        for v_neighbor in adj[u]:
            if v_neighbor not in visited:
                visited.add(v_neighbor)
                q.append(v_neighbor)
    
    return len(visited) == n_vertices

def find_n_for_cycle_graph(limit):
    """Finds all n up to a limit for which AG(Z_n) is a cycle graph."""
    solutions = []
    for n in range(2, limit + 1):
        graph = build_associate_graph(n)
        if is_cycle_graph(graph):
            solutions.append(n)
    return solutions

# --- Main Execution ---
# Set a reasonable limit to test for n.
# Our mathematical proof shows no solutions exist, this code verifies it computationally.
SEARCH_LIMIT = 50
solutions = find_n_for_cycle_graph(SEARCH_LIMIT)

# Format the output as requested by the user prompt.
if not solutions:
    # As expected, the list of solutions is empty.
    print("n \\in\\{\\}")
else:
    # This case should not be reached.
    n_str = ", ".join(map(str, solutions))
    print(f"n \\in\\{{{n_str}\\}}")