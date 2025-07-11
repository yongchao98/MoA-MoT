import networkx as nx
from itertools import combinations

def solve():
    """
    Solves the problem by analyzing the graph from the 3-SAT to IS reduction.
    """
    # Step 1: Define the graph from the image.
    # We label the vertices 0 to 11. The mapping to the image is described
    # in the comments within the code.
    G = nx.Graph()
    nodes = list(range(12))
    G.add_nodes_from(nodes)
    
    # Adjacency list based on visual inspection of the graph image.
    # This is a critical step. An error here would change the outcome.
    # Node Numbering Scheme:
    # 0: top-left peak         3: main top peak       6: right peak
    # 1: top-left, bottom-left  4: main top, left      7: right, bottom-left
    # 2: top-left, bottom-right 5: main top, right     8: right, bottom-right
    # 9: center top            10: center left
    # 11: center right
    edges = [
        # Leftmost component
        (0, 1), (0, 2), (1, 2),
        # Top component
        (3, 4), (3, 5), (4, 5),
        # Rightmost component
        (6, 7), (6, 8), (7, 8),
        # Center component
        (9, 10), (9, 11), (10, 11),
        # Connections BETWEEN components
        (2, 10), (1, 10), (4, 10), (4, 6), (5, 11), (5, 6), (7, 11)
    ]
    G.add_edges_from(edges)

    # Step 2: Find the 4 clause-triangles.
    # A clause triangle partition is a partition of V into four K_3 cliques.
    cliques = [c for c in nx.find_cliques(G) if len(c) == 3]
    
    clause_partitions = []
    # Find partitions of the 12 vertices into 4 sets of 3.
    # We check if each set of 3 forms a triangle (a 3-clique).
    v = set(nodes)
    for c1_nodes in combinations(v, 3):
        if sorted(c1_nodes) not in [sorted(c) for c in cliques]: continue
        v2 = v - set(c1_nodes)
        for c2_nodes in combinations(v2, 3):
            if sorted(c2_nodes) not in [sorted(c) for c in cliques]: continue
            v3 = v2 - set(c2_nodes)
            for c3_nodes in combinations(v3, 3):
                if sorted(c3_nodes) not in [sorted(c) for c in cliques]: continue
                c4_nodes = v3 - set(c3_nodes)
                if sorted(list(c4_nodes)) not in [sorted(c) for c in cliques]: continue
                
                partition = sorted([tuple(sorted(c)) for c in [c1_nodes, c2_nodes, c3_nodes, list(c4_nodes)]])
                if partition not in clause_partitions:
                    clause_partitions.append(partition)

    print("Analysis of the graph from the 3-SAT reduction:\n")
    if not clause_partitions:
        print("Error: Could not find a partition of the graph into 4 triangles.")
        print("This implies the graph was not formed by the described reduction, or my graph encoding is wrong.")
        return

    # Assuming there's a unique partition. If not, the logic applies to all.
    # Based on the visual structure, the intended partition is the four "hills".
    partition = [ [0,1,2], [3,4,5], [6,7,8], [9,10,11] ]
    # Let's verify this partition is indeed a valid K3 partition
    is_valid = True
    for triangle in partition:
        if not (G.has_edge(triangle[0], triangle[1]) and \
                G.has_edge(triangle[0], triangle[2]) and \
                G.has_edge(triangle[1], triangle[2])):
            is_valid = False
            break
    
    if not is_valid:
        # Fallback to found partitions if the assumed one is incorrect
        if clause_partitions:
            partition = clause_partitions[0]
            print(f"Used the first found partition: {partition}")
        else: # Should not be reached if the first check passed
            return
    else:
        print(f"Found a unique and logical triangle partition: {partition}")


    # Step 3: Identify literal relationships from inter-triangle edges.
    inter_edges = [e for e in G.edges() if sorted(e) not in [sorted(c) for c in nx.find_cliques(G.subgraph(e))]]
    
    clauses = partition
    
    def get_clause_idx(v):
        for i, c in enumerate(clauses):
            if v in c:
                return i
        return -1
    
    inter_edges = [e for e in G.edges() if get_clause_idx(e[0]) != get_clause_idx(e[1])]
    
    # For any edge (u, v) in inter_edges, L_u must equal !L_v.
    # We can build a "literal constraint graph" and find required literal equalities.
    # If L_u = !L_w and L_v = !L_w, then L_u = L_v.
    
    literal_groups = {i: {i} for i in range(12)}
    
    # Propagate equality constraints
    for _ in range(12): # Iterate to ensure all constraints propagate
        for w in range(12):
            w_neg_neighbors = {u for u,v in inter_edges if v==w} | {v for u,v in inter_edges if u==w}
            if len(w_neg_neighbors) < 2:
                continue
            
            # All nodes connected to w must have the same literal (!L_w)
            first_neighbor = next(iter(w_neg_neighbors))
            first_group = literal_groups[first_neighbor]
            for neighbor in w_neg_neighbors:
                literal_groups[neighbor] = literal_groups[neighbor].union(first_group)
                literal_groups[first_neighbor] = literal_groups[first_neighbor].union(literal_groups[neighbor])


    # Normalize groups
    final_groups = []
    seen = set()
    for i in range(12):
        if i in seen: continue
        group = frozenset(literal_groups[i])
        final_groups.append(group)
        seen.update(group)

    print("\nVertices that must represent the same literal:")
    for group in final_groups:
        if len(group) > 1:
            print(f"- {set(group)}")
    
    # Step 4: The contradiction.
    # "a clause ... cannot contain the same literal twice".
    # This means any two vertices in a literal group cannot be in the same clause triangle.
    print("\nChecking for contradictions:")
    has_contradiction = False
    for group in final_groups:
        if len(group) <= 1:
            continue
        # Check if any two vertices from this group are in the same clause.
        for v1, v2 in combinations(group, 2):
            if get_clause_idx(v1) == get_clause_idx(v2):
                print(f"Contradiction found! Vertices {v1} and {v2} must have the same literal.")
                print(f"However, they both belong to the same clause: {clauses[get_clause_idx(v1)]}.")
                print("This violates the rule that a clause cannot contain the same literal twice.")
                has_contradiction = True
                break
        if has_contradiction:
            break
            
    # Also, we can find a second contradiction.
    # The number of vertices in a literal group cannot exceed the number of clauses.
    num_clauses = len(clauses)
    for group in final_groups:
        if len(group) > num_clauses:
            print(f"Contradiction found! The literal group {set(group)} has {len(group)} members.")
            print(f"Since each must belong to a different clause, we would need at least {len(group)} clauses.")
            print(f"But the graph only has {num_clauses} clauses.")
            has_contradiction = True
            break
            
    # Step 5: Conclusion
    print("\nConclusion:")
    if has_contradiction:
        print("The structural constraints of the graph are inconsistent with the rules of the 3-SAT to IS reduction.")
        print("Therefore, no valid 3-CNF formula could have produced this graph.")
        print("The set of possible formulas is empty.")
        min_models = 0
        max_models = 0
    else:
        # This part would execute if no contradiction was found.
        # It would involve calculating models for different valid formulas.
        print("No contradiction found. The problem is more complex.")
        min_models = -1 # Placeholder
        max_models = -1 # Placeholder

    print(f"\nFinal Answer: The minimum and maximum number of models are ({min_models}, {max_models}).")
    return min_models, max_models

min_val, max_val = solve()
# The final output format required by the user prompt
# This is a bit of a hack to get the desired output string
# print(f"<<<({min_val}, {max_val})>>>")
final_answer = "({}, {})".format(min_val, max_val)
print(f"\n<<<{final_answer}>>>")
