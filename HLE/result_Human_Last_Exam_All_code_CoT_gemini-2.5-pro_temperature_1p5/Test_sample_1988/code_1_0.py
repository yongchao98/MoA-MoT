import networkx as nx
import time

def solve():
    """
    This function calculates the number of subgraphs with HoG ID 50698
    (the Kneser graph K(8,2)) contained in the Gosset graph.
    """
    print("Step 1: Constructing the graphs...")
    # The subgraph is the Kneser graph K(8,2), as per the description.
    # Vertices are the C(8,2)=28 edges of K8.
    # Adjacency if edges are disjoint.
    subgraph = nx.kneser_graph(8, 2)
    
    # The target graph is the Gosset graph.
    target_graph = nx.gosset_graph()
    
    print(f"Target Graph (Gosset): {target_graph.number_of_vertices()} vertices, {target_graph.number_of_edges()} edges.")
    print(f"Subgraph (Kneser K(8,2)): {subgraph.number_of_vertices()} vertices, {subgraph.number_of_edges()} edges.")
    print("-" * 30)

    # --- Find the number of automorphisms of the subgraph ---
    # An automorphism is an isomorphism of a graph to itself.
    print("Step 2: Calculating the number of automorphisms for the subgraph K(8,2)...")
    print("(This may take a moment)")
    start_time_aut = time.time()
    
    # We use GraphMatcher to find all isomorphisms from the subgraph to itself.
    aut_matcher = nx.isomorphism.GraphMatcher(subgraph, subgraph)
    num_automorphisms = 0
    # The number of isomorphisms is the size of the automorphism group.
    for _ in aut_matcher.isomorphisms_iter():
        num_automorphisms += 1
        
    end_time_aut = time.time()
    print(f"Found {num_automorphisms} automorphisms in {end_time_aut - start_time_aut:.2f} seconds.")
    print("-" * 30)

    # --- Find the number of monomorphisms (isomorphisms to induced subgraphs) ---
    print("Step 3: Calculating the total number of induced subgraphs matching K(8,2) in the Gosset graph...")
    print("(This is computationally intensive and may take several minutes)")
    start_time_mono = time.time()
    
    # We use GraphMatcher to find all monomorphisms (isomorphisms to induced subgraphs).
    mono_matcher = nx.isomorphism.GraphMatcher(target_graph, subgraph)
    num_monomorphisms = 0
    # We must check for induced subgraphs specifically.
    for _ in mono_matcher.subgraph_isomorphisms_iter(induced=True):
        num_monomorphisms += 1
        
    end_time_mono = time.time()
    print(f"Found {num_monomorphisms} total monomorphisms (isomorphisms to subgraphs) in {end_time_mono - start_time_mono:.2f} seconds.")
    print("-" * 30)

    # --- Calculate the final result ---
    print("Step 4: Calculating the number of unique subgraphs...")
    # The number of unique subgraphs is the number of monomorphisms divided by the number of automorphisms.
    # Use integer division as the result must be an integer.
    if num_automorphisms > 0:
        num_unique_subgraphs = num_monomorphisms // num_automorphisms
    else:
        num_unique_subgraphs = 0
        
    print("Final Calculation:")
    print(f"Number of unique subgraphs = {num_monomorphisms} (monomorphisms) / {num_automorphisms} (automorphisms) = {num_unique_subgraphs}")

solve()
<<<28>>>