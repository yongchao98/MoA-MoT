import math
from collections import deque

def get_units(n):
    """Returns the set of units in Z_n."""
    if n == 1:
        return {0}
    return {i for i in range(1, n) if math.gcd(i, n) == 1}

def build_associate_graph(n):
    """Builds the associate graph AG(Z_n) and returns its adjacency list."""
    if n <= 1:
        return {}
    
    vertices = list(range(1, n))
    units = get_units(n)
    adj = {v: [] for v in vertices}

    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            a = vertices[i]
            b = vertices[j]
            # Check if a and b are associates
            # a = b*u (mod n) for some unit u.
            # This is equivalent to checking if a is in the orbit of b.
            # A simple way is to check if b*u == a for any unit u.
            for u in units:
                if (b * u) % n == a:
                    adj[a].append(b)
                    adj[b].append(a)
                    break # Found associate, no need to check other units
    return adj

def is_cycle(graph, n):
    """Checks if a graph is a cycle graph."""
    num_vertices = n - 1
    if num_vertices < 3:
        return False

    # 1. Every vertex must have degree 2
    for v in graph:
        if len(graph[v]) != 2:
            return False

    # 2. The graph must be connected
    q = deque()
    q.append(next(iter(graph))) # Start BFS from an arbitrary vertex
    visited = {next(iter(graph))}
    count = 0
    while q:
        v = q.popleft()
        count += 1
        for neighbor in graph[v]:
            if neighbor not in visited:
                visited.add(neighbor)
                q.append(neighbor)
                
    if count != num_vertices:
        return False

    return True

def find_n_for_ring_graph(limit=50):
    """Finds all n <= limit for which AG(Z_n) is a ring graph."""
    solutions = []
    # A cycle must have >= 3 vertices, so n-1 >= 3 => n >= 4
    for n in range(4, limit + 1):
        graph = build_associate_graph(n)
        if is_cycle(graph, n):
            solutions.append(n)
    return solutions

if __name__ == "__main__":
    solutions = find_n_for_ring_graph(100)
    
    # As per the prompt, print the final answer in the requested format
    result_str = ", ".join(map(str, solutions))
    print(f"n in{{ {result_str} }}")
