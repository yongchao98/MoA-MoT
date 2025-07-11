import math
from collections import deque

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    while b:
        a, b = b, a % b
    return a

def is_cycle(adj, num_vertices):
    """
    Checks if a simple graph is a cycle.
    A cycle must be connected, have all vertices with degree 2,
    and must have at least 3 vertices.
    """
    # A cycle must have at least 3 vertices.
    if num_vertices < 3:
        return False
    
    all_vertices = set(adj.keys())
    
    # Check if any vertex has a degree other than 2.
    for i in range(1, num_vertices + 1):
         if adj.get(i, []) != 2:
            return False

    # Check for connectivity using Breadth-First Search (BFS).
    # If the graph is not connected, the BFS will not visit all vertices.
    q = deque()
    start_node = 1
    q.append(start_node)
    
    visited = {start_node}
    
    while q:
        node = q.popleft()
        for neighbor in adj.get(node, []):
            if neighbor not in visited:
                visited.add(neighbor)
                q.append(neighbor)
                
    return len(visited) == num_vertices

def find_associate_ring_graph_cycles(limit):
    """
    Finds all integers n up to a given limit for which the associate
    ring graph AG(Z_n) is a cycle.
    """
    cycle_ns = []
    # A cycle must have at least 3 vertices, so n-1 >= 3, which means n >= 4.
    for n in range(4, limit + 1):
        vertices = set(range(1, n))
        num_vertices = n - 1
        
        # Find the set of units in Z_n
        units = {i for i in range(1, n) if gcd(i, n) == 1}
        
        adj = {v: [] for v in vertices}
        adj_degrees = {v: 0 for v in vertices}

        unclassified_vertices = set(vertices)
        
        while unclassified_vertices:
            v_seed = unclassified_vertices.pop()
            
            # Calculate the associate class for the seed vertex
            current_class = set()
            for u in units:
                associate = (v_seed * u) % n
                if associate in vertices: # Make sure not to include the 0 element
                    current_class.add(associate)

            # Build a clique on the vertices of the current class
            class_list = sorted(list(current_class))
            for i in range(len(class_list)):
                for j in range(i + 1, len(class_list)):
                    v1, v2 = class_list[i], class_list[j]
                    adj_degrees[v1] += 1
                    adj_degrees[v2] += 1

            # Remove all members of the class from the unclassified set
            unclassified_vertices -= current_class

        # Check if the graph is a cycle.
        # First, a quick check on degrees. All must be 2.
        is_2_regular = all(degree == 2 for degree in adj_degrees.values())
        
        if is_2_regular:
            # Reconstruct adjacency list for connectivity check
            adj_list = {v: [] for v in vertices}
            # This part can be combined with the class processing loop,
            # but is separated for clarity.
            unclassified_vertices_2 = set(vertices)
            while unclassified_vertices_2:
                v_seed = unclassified_vertices_2.pop()
                current_class = set()
                for u in units:
                    associate = (v_seed * u) % n
                    if associate in vertices:
                        current_class.add(associate)
                class_list = sorted(list(current_class))
                for i in range(len(class_list)):
                    for j in range(i + 1, len(class_list)):
                        v1, v2 = class_list[i], class_list[j]
                        adj_list[v1].append(v2)
                        adj_list[v2].append(v1)
                unclassified_vertices_2 -= current_class

            # A 2-regular graph is a disjoint union of cycles.
            # To be a single cycle, it must be connected.
            if is_cycle(adj_list, num_vertices):
                 cycle_ns.append(n)
            
    return cycle_ns

# Set a search limit and find the solutions.
# Based on the mathematical analysis, this list is expected to be empty.
limit = 100
solutions = find_associate_ring_graph_cycles(limit)

# Format the final output as requested.
solution_str = ", ".join(map(str, solutions))
print(f"n \\in{{{solution_str}}}")
>>>n \in{}