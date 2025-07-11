import networkx as nx
from itertools import combinations

def count_slices(graph_name, G):
    """
    Calculates and explains the number of slices N(G) for a given cubic graph G.
    """
    nodes = list(G.nodes())
    edges = [tuple(sorted(e)) for e in G.edges()]
    num_edges = len(edges)
    
    adj = {node: set() for node in nodes}
    for u, v in edges:
        adj[u].add(tuple(sorted((u,v))))
        adj[v].add(tuple(sorted((u,v))))
        
    slice_count = 0
    
    partition_details = {}

    # Iterate through partitions of edges. A partition is defined by one set.
    # To avoid counting pairs (A, B) and (B, A) twice, we iterate sizes up to num_edges/2.
    for i in range(1, num_edges // 2 + 1):
        for class1_edges_tuple in combinations(edges, i):
            class1_edges = set(class1_edges_tuple)
            
            is_a_slice = True
            for node in nodes:
                c1_incident = False
                c2_incident = False
                for edge in adj[node]:
                    if edge in class1_edges:
                        c1_incident = True
                    else: # edge is in class2
                        c2_incident = True
                
                if not (c1_incident and c2_incident):
                    is_a_slice = False
                    break
            
            if is_a_slice:
                partition_size_pair = tuple(sorted((len(class1_edges), num_edges - len(class1_edges))))
                if partition_size_pair not in partition_details:
                    partition_details[partition_size_pair] = 0
                partition_details[partition_size_pair] += 1
                slice_count +=1
                
    print(f"For {graph_name} (m={len(nodes)}, |E|={num_edges}):")
    # Present the breakdown of slice calculations
    equation_parts = []
    for sizes, count in sorted(partition_details.items()):
        equation_parts.append(f"{count}")
        print(f"   Found {count} slices from partitions of sizes {sizes}.")
    
    if len(equation_parts) > 1:
        equation_str = " + ".join(equation_parts)
        print(f"   Total N({graph_name}) = {equation_str} = {slice_count}.")
    else:
        print(f"   Total N({graph_name}) = {slice_count}.")
        
    return slice_count


# --- Main logic to determine M(n) values ---

# M(0)
print("Step 1: Determine M(0)")
print("N(G) must be a multiple of 0, which means N(G)=0.")
print("It is a long-standing mathematical conjecture that every simple cubic graph has a slice, so N(G) > 0 for all G.")
print("Assuming this conjecture, no such graph exists.")
m0_answer = "none"
print("Result: M(0) is none.")
print("-" * 20)

# M(3)
print("Step 2: Determine M(3)")
print("Searching for the smallest m where N(G) is a multiple of 3.")
print("The smallest cubic graph has m=4 vertices, which is K_4 (the complete graph on 4 vertices).")
k4 = nx.complete_graph(4)
n_k4 = count_slices("K_4", k4)
if n_k4 % 3 == 0:
    m3_answer = 4
    print(f"N(K_4)={n_k4} is a multiple of 3. Since m=4 is the smallest possible, M(3) = 4.")
else:
    # This branch is not expected based on manual calculation.
    m3_answer = "unknown"
    print(f"N(K_4)={n_k4} is not a multiple of 3. Would need to check m=6.")
print("-" * 20)


# M(5)
print("Step 3: Determine M(5)")
print("Searching for the smallest m where N(G) is a multiple of 5.")
print(f"For m=4, N(K_4) = {n_k4}, which is not a multiple of 5.")
print("The next smallest size for a cubic graph is m=6.")
print("Checking the prism graph, one of the two cubic graphs on 6 vertices.")
prism_graph = nx.Graph()
# Prism graph vertices are 1-indexed for clarity matching drawings
prism_graph.add_edges_from([(1,2),(2,3),(3,1), (4,5),(5,6),(6,4), (1,4),(2,5),(3,6)])
n_prism = count_slices("Prism", prism_graph)
m5_answer = "unknown"
if n_prism % 5 == 0:
    m5_answer = 6
    print(f"N(Prism)={n_prism} is a multiple of 5. As m=4 did not work, the smallest m is 6. M(5) = 6.")
else:
    # Check the other m=6 graph if needed
    k33 = nx.complete_bipartite_graph(3,3)
    k33_nodes = {i: i+1 for i in range(6)} # relabel to 1-6 for consistency
    k33 = nx.relabel_nodes(k33,k33_nodes)
    n_k33 = count_slices("K_3,3", k33)
    if n_k33 % 5 == 0:
      m5_answer = 6
      print(f"N(K_3,3)={n_k33} is a multiple of 5. As m=4 did not work, the smallest m is 6. M(5) = 6.")

print("-" * 20)
# Final Answer Formatting
final_answer_string = f"{m0_answer},{m3_answer},{m5_answer}"
print(f"Final determined values: M(0)={m0_answer}, M(3)={m3_answer}, M(5)={m5_answer}")
print(f"Final answer in required format: {final_answer_string}")
print(f'<<<{final_answer_string}>>>')
