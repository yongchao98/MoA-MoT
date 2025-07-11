import collections

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_graph_a_cycle(n, edges, num_vertices):
    """
    Checks if a graph is a cycle.
    A graph is a cycle if it is connected and every vertex has a degree of 2.
    """
    if num_vertices < 3:
        return False

    degrees = collections.defaultdict(int)
    adj = collections.defaultdict(list)
    for u, v in edges:
        degrees[u] += 1
        degrees[v] += 1
        adj[u].append(v)
        adj[v].append(u)

    # Check if all vertices have degree 2
    for i in range(1, n):
        if degrees[i] != 2:
            return False

    # Check for connectivity
    # A 2-regular graph is connected if and only if it doesn't have multiple components.
    # Since all vertices have degree 2, the graph is a disjoint union of cycles.
    # For it to be a single cycle, it must be connected.
    # We can do a traversal (like BFS or DFS) from an arbitrary node to see if we visit all nodes.
    q = collections.deque([1])
    visited = {1}
    count = 0
    while q:
        node = q.popleft()
        count += 1
        for neighbor in adj[node]:
            if neighbor not in visited:
                visited.add(neighbor)
                q.append(neighbor)
    
    return count == num_vertices

def find_n_for_ring_graph(search_limit):
    """
    Finds values of n for which AG(Z_n) is a ring graph (cycle).
    """
    solutions = []
    
    for n in range(4, search_limit + 1):
        num_vertices = n - 1
        vertices = list(range(1, n))
        
        # Find all units in Z_n
        units = [u for u in vertices if gcd(u, n) == 1]
        
        # Build the graph
        edges = []
        for i in range(num_vertices):
            for j in range(i + 1, num_vertices):
                v1 = vertices[i]
                v2 = vertices[j]
                
                # Check if v1 and v2 are associates
                # v1 = v2 * u (mod n) for some unit u
                for u in units:
                    if (v2 * u) % n == v1:
                        edges.append((v1, v2))
                        break # Found an association, no need to check other units
        
        if is_graph_a_cycle(n, edges, num_vertices):
            solutions.append(n)
            
    return solutions

# As derived in the logical plan, no such integers n exist.
# The code will confirm this by finding no solutions.
# We run the search up to a reasonable limit.
found_solutions = find_n_for_ring_graph(100)

# The question asks to write all values of integers in a sequence.
# Since the set of solutions is empty, we represent it as an empty set.
if not found_solutions:
    print("n \u2208 {\u007D") # n ∈ {}
else:
    # This part of the code is not expected to be reached.
    sequence = ", ".join(map(str, found_solutions))
    print(f"n \u2208 {{ {sequence} }}")
