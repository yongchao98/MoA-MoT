import itertools

def count_homomorphisms(query_graph, target_graph):
    """
    Counts the number of homomorphisms from a query_graph to a target_graph.
    Graphs are represented as adjacency lists (dictionaries).
    This is a brute-force implementation and is slow for large graphs.
    """
    q_nodes = sorted(query_graph.keys())
    t_nodes = sorted(target_graph.keys())
    
    if not q_nodes:
        return 1

    count = 0
    # Iterate through all possible functions from V(query) to V(target)
    for mapping_tuple in itertools.product(t_nodes, repeat=len(q_nodes)):
        # f is the current mapping function
        f = {q_nodes[i]: mapping_tuple[i] for i in range(len(q_nodes))}
        
        is_homomorphism = True
        # Check the homomorphism condition for every edge in the query graph
        for u, neighbors in query_graph.items():
            for v in neighbors:
                # To handle directed vs undirected, we check one direction is enough
                # if the query graph representation is symmetric.
                if u > v:
                    continue
                # Map the edge (u, v) to (f(u), f(v))
                mapped_u = f[u]
                mapped_v = f[v]
                
                # Check if (f(u), f(v)) is an edge in the target graph
                if mapped_v not in target_graph.get(mapped_u, []):
                    is_homomorphism = False
                    break
            if not is_homomorphism:
                break
        
        if is_homomorphism:
            count += 1
            
    return count

# Define a target graph G (a triangle with a tail)
G = {
    0: [1, 2],
    1: [0, 2],
    2: [0, 1, 3],
    3: [2]
}

# Define an acyclic query represented by a forest Q.
# Q is the disjoint union of a path of length 1 (T1) and a single vertex (T2).
T1 = {'a': ['b'], 'b': ['a']} # A single edge (K_2)
T2 = {'c': []}               # A single vertex (K_1)
Q = {'a': ['b'], 'b': ['a'], 'c': []} # The forest T1 U T2

# Calculate the number of answers (homomorphisms)
hom_T1_G = count_homomorphisms(T1, G)
hom_T2_G = count_homomorphisms(T2, G)
hom_Q_G = count_homomorphisms(Q, G)

product_of_tree_homs = hom_T1_G * hom_T2_G

print("This code demonstrates that for a forest Q = T1 U T2, hom(Q, G) = hom(T1, G) * hom(T2, G).")
print("-" * 20)
print(f"Target Graph G has {len(G)} vertices and {sum(len(v) for v in G.values()) // 2} edges.")
print(f"Query is a forest Q with {len(Q)} vertices.")
print(f"It consists of tree T1 (edge) and T2 (vertex).")
print("-" * 20)
print("Number of homomorphisms from T1 to G:")
print(f"hom(T1, G) = {hom_T1_G}")
print("\nNumber of homomorphisms from T2 to G:")
print(f"hom(T2, G) = {hom_T2_G}")
print("\nProduct of the two homomorphism counts:")
print(f"hom(T1, G) * hom(T2, G) = {hom_T1_G} * {hom_T2_G} = {product_of_tree_homs}")
print("\nNumber of homomorphisms from the forest Q to G (calculated directly):")
print(f"hom(Q, G) = {hom_Q_G}")
print("-" * 20)

if hom_Q_G == product_of_tree_homs:
    print(f"As shown, {hom_Q_G} = {product_of_tree_homs}. The number of answers for the forest query equals the product of answers for its component trees.")
    print("This confirms the logic used in the proof.")
else:
    print("There is an error in the logic or calculation.")
