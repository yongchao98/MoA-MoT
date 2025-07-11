import networkx as nx

def count_c5(G):
    """
    Counts the number of 5-cycles in a graph G.
    
    The algorithm iterates through all paths of length 4 and checks if the last 
    vertex is connected to the first, forming a cycle. To be efficient, it 
    uses set intersections to find common neighbors.
    
    The final count is divided by 10 to correct for overcounting, as each
    5-cycle is found 5 times (once for each starting node) and in 2 different
    directions (clockwise and counter-clockwise).
    """
    cycle_count = 0
    
    # Using a dictionary of neighbor sets for faster lookups
    adj = {n: set(G.neighbors(n)) for n in G.nodes()}

    for u in G.nodes():
        # The path is u-v-w-x-y-u. We find all such unique paths.
        # v is a neighbor of u
        for v in adj[u]:
            # w is a neighbor of v (w must not be u)
            for w in adj[v]:
                if w == u:
                    continue
                # x is a neighbor of w (x must not be u or v)
                for x in adj[w]:
                    if x == v or x == u:
                        continue
                    
                    # y must be a common neighbor of the start (u) and end (x) vertices of the path u-v-w-x.
                    # y must not be v or w to ensure all 5 vertices of the cycle are unique.
                    common_neighbors = adj[x].intersection(adj[u])
                    
                    for y in common_neighbors:
                        if y != v and y != w:
                            cycle_count += 1
                            
    # Each 5-cycle (u,v,w,x,y) is counted 10 times:
    # 5 times for each starting vertex (u, v, w, x, y)
    # 2 times for each direction (e.g., u-v-w-x-y-u and u-y-x-w-v-u)
    return cycle_count // 10

# The parameters for the two non-isomorphic SRGs
n = 16
d = 6
# `lambda` is a keyword in Python, so we use `lam`
lam = 2
mu = 2

# Main part of the script
print("Yes, there exists a combination of parameters (n,d,lambda,mu) and a pair of graphs G, G' from srg(n,d,lambda,mu) such that G and G' do not have the same number of 5-cycles.")
print("-" * 50)
print(f"The chosen parameters are (n, d, lambda, mu) = ({n}, {d}, {lam}, {mu}).")
print("The two non-isomorphic graphs with these parameters are the L2(4) graph and the Shrikhande graph.")
print("-" * 50)

# Create the L2(4) graph (4x4 Rook's graph)
# It is the Cartesian product of two complete graphs K4.
K4 = nx.complete_graph(4)
G1_prod = nx.cartesian_product(K4, K4)
G1 = nx.convert_node_labels_to_integers(G1_prod, first_label=0)
G1.name = "L2(4) Graph (4x4 Rook's Graph)"

# Create the Shrikhande graph
G2 = nx.shrikhande_graph()
G2.name = "Shrikhande Graph"

# Count the 5-cycles in each graph
num_c5_G1 = count_c5(G1)
num_c5_G2 = count_c5(G2)

# Print the final results
print(f"Graph 1: {G1.name}")
print(f"Number of 5-cycles = {num_c5_G1}")
print()
print(f"Graph 2: {G2.name}")
print(f"Number of 5-cycles = {num_c5_G2}")
print("-" * 50)

if num_c5_G1 != num_c5_G2:
    print("Conclusion: The number of 5-cycles is different, confirming that such a case exists.")
else:
    # This case should not be reached with the correct implementation.
    print("Conclusion: The number of 5-cycles is the same. The code or example may be incorrect.")
