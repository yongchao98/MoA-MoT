import math
from collections import defaultdict, deque

def get_proper_divisors(n):
    """Returns a list of proper divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            if i < n:
                divs.add(i)
            if n // i < n:
                divs.add(n // i)
    return sorted(list(divs))

def euler_phi(n):
    """Computes Euler's totient function."""
    if n == 0: return 0
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return int(result)

def is_cycle(adj, num_vertices):
    """Checks if a graph is a cycle."""
    if num_vertices < 3:
        return False
    
    # Every vertex must have degree 2
    for v in range(1, num_vertices + 1):
        if v not in adj or len(adj[v]) != 2:
            return False

    # Check for connectivity and single component
    q = deque()
    visited = set()
    
    start_node = 1
    q.append(start_node)
    visited.add(start_node)
    
    count = 0
    while q:
        u = q.popleft()
        count += 1
        for v in adj[u]:
            if v not in visited:
                visited.add(v)
                q.append(v)
    
    return count == num_vertices

def find_ring_graph_n():
    """
    Finds all integers n for which AG(Z_n) is a ring graph.
    """
    solutions = []
    # A cycle needs at least 3 vertices, so n-1 >= 3 => n >= 4.
    for n in range(4, 151): # Check a reasonable range for n.
        num_vertices = n - 1
        adj = defaultdict(list)
        vertices = list(range(1, n))
        
        divisors = get_proper_divisors(n)
        
        # Construct the graph based on the disjoint cliques model
        for d in divisors:
            # All elements 'x' with gcd(x,n)=d form a clique.
            # The size of this clique is phi(n/d).
            clique_nodes = [v for v in vertices if math.gcd(v, n) == d]
            if len(clique_nodes) > 1:
                for i in range(len(clique_nodes)):
                    for j in range(i + 1, len(clique_nodes)):
                        u, v = clique_nodes[i], clique_nodes[j]
                        adj[u].append(v)
                        adj[v].append(u)

        if is_cycle(adj, num_vertices):
            solutions.append(n)

    # Based on the mathematical derivation, this list will be empty.
    # The output format follows the problem's request.
    result_str = ",".join(map(str, solutions))
    print(f"n \u2208 {{{result_str}}}")

# Execute the search and print the result
find_ring_graph_n()