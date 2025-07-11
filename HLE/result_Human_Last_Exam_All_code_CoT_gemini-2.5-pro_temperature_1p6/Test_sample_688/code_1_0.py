import sys
import itertools
from math import comb, factorial
from fractions import Fraction

def is_connected(n_vertices_subgraph, adj, nodes_to_check, start_node, removed_node=None):
    """Check if the graph represented by adj is connected using BFS."""
    if n_vertices_subgraph == 0:
        return True
    
    q = [start_node]
    visited = {start_node}
    
    head = 0
    while head < len(q):
        u = q[head]
        head += 1
        
        # In this implementation, adj is the adjacency list for the full graph K_n.
        # We only consider neighbors that are part of the current subgraph.
        for v in adj[u]:
            if v != removed_node and v not in visited:
                visited.add(v)
                q.append(v)
                
    return len(visited) == n_vertices_subgraph

def is_biconnected(n_vertices, edges):
    """
    Check if a graph is biconnected.
    A graph is biconnected if it is connected and has no articulation points.
    We assume the graph is on n_vertices and some nodes might be isolated.
    """
    if n_vertices < 3:
        return False

    adj = [[] for _ in range(n_vertices)]
    nodes_with_edges = set()
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        nodes_with_edges.add(u)
        nodes_with_edges.add(v)

    # All vertices must be part of the graph for it to be biconnected on n vertices
    if len(nodes_with_edges) != n_vertices:
        return False

    # 1. Check for global connectivity
    if not is_connected(n_vertices, adj, list(range(n_vertices)), 0):
        return False
        
    # 2. Check for articulation points
    for i in range(n_vertices):
        # Temporarily remove vertex i
        remaining_nodes = [node for node in range(n_vertices) if node != i]
        if not remaining_nodes:
             continue
        
        # Check connectivity of the remaining graph
        if not is_connected(n_vertices - 1, adj, remaining_nodes, remaining_nodes[0], removed_node=i):
            return False # Found an articulation point
            
    return True

def calculate_cn(n):
    """
    Calculates the prefactor c_n for the fully f-connected Ree-Hoover diagram.
    """
    if n < 2:
        raise ValueError("n must be at least 2.")
    
    if n == 2:
        # For n=2, B_2 = -1/2 * integral(f_12), and Lambda_2 = integral(f_12).
        # Thus, c_2 = -1/2.
        return Fraction(-1, 2)

    vertices = range(n)
    possible_edges = list(itertools.combinations(vertices, 2))
    num_total_edges = len(possible_edges)

    w_Kn = 0
    
    # Iterate through all possible subgraphs H of K_n by number of edges
    for num_edges in range(1, num_total_edges + 1):
        for edges in itertools.combinations(possible_edges, num_edges):
            if is_biconnected(n, list(edges)):
                term = (-1)**(num_total_edges - num_edges)
                w_Kn += term

    cn_numerator = -(n - 1) * w_Kn
    cn_denominator = factorial(n)
    
    return Fraction(cn_numerator, cn_denominator)

if __name__ == '__main__':
    print("Calculating the prefactor c_n for n = 3, 4, 5.")
    print("The prefactor appears in the virial expansion as:")
    print("B_n = c_n * Lambda_n + B_n^ex")
    print("-" * 30)

    for n_val in range(3, 6):
        c_n = calculate_cn(n_val)
        print(f"For n = {n_val}:")
        # To match the requested format "output each number in the final equation"
        print(f"B_{n_val} = ({c_n.numerator}/{c_n.denominator}) * Lambda_{n_val} + B_{n_val}^ex")

    # Example calculation for a user-provided n via command line argument
    if len(sys.argv) > 1:
        try:
            n_input = int(sys.argv[1])
            if n_input < 2:
                print(f"\nError: n must be an integer >= 2.")
            elif n_input > 5:
                print(f"\nCalculating for user input n = {n_input} (this might take a while)...")
                c_n = calculate_cn(n_input)
                print(f"B_{n_input} = ({c_n.numerator}/{c_n.denominator}) * Lambda_{n_input} + B_{n_input}^ex")
        except ValueError:
            print(f"\nCould not parse command-line argument '{sys.argv[1]}' as an integer.")
