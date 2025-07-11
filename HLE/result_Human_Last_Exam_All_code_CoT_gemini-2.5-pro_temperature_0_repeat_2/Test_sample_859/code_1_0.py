import math

def solve():
    """
    This function calculates the minimal number of new edges to make G' 2-edge-connected.
    The reasoning is based on a worst-case construction.
    
    Let k be the number of connected components of G', and l be the total number of leaves in its condensation forest.
    The number of edges to add is m = (k-1) + ceil((l - 2*(k-1))/2).
    
    We can construct a graph G that satisfies the given conditions, for which the resulting G'
    has k=3 components, and each component is a 2-edge-connected block (i.e., a leaf in the forest).
    So for this construction, k=3 and l=3.
    
    Let's calculate m for this case.
    """
    
    k = 3  # Number of components in the constructed G'
    l = 3  # Number of leaves in the condensation forest of the constructed G'
    
    # We need k-1 edges to connect the components
    edges_to_connect = k - 1
    
    # After connecting, we have one component. The number of leaves in the new condensation tree is l' = l - 2*(k-1)
    # This is because each of the k-1 new edges connects two leaves, reducing the total leaf count by 2.
    leaves_after_connecting = l - 2 * edges_to_connect
    
    # We need ceil(l'/2) edges to make the resulting single tree 2-edge-connected.
    # Note: if leaves_after_connecting is negative, it means we have created cycles and have fewer than 2 leaves.
    # The ceil function handles this correctly, resulting in 0 or fewer additional edges.
    if leaves_after_connecting < 0:
        edges_for_2_connectivity = math.ceil(leaves_after_connecting / 2.0)
        # In graph theory terms, adding an edge between two leaves of a tree reduces the leaf count by 2.
        # If we add k-1 edges to connect k trees, we can reduce the leaf count by at most 2(k-1).
        # If l < 2(k-1), we can form a single component with cycles and no leaves, which is 2-edge-connected.
        # The formula can result in a negative value inside ceil, which means we have "over-connected" the graph.
        # A more standard interpretation is that we need max(0, ceil(l'/2)) edges.
        edges_for_2_connectivity = max(0, math.ceil(leaves_after_connecting / 2.0))

    else:
        edges_for_2_connectivity = math.ceil(leaves_after_connecting / 2.0)

    total_edges = edges_to_connect + edges_for_2_connectivity
    
    print("Worst-case construction:")
    print(f"Number of components (k) = {k}")
    print(f"Number of leaves (l) = {l}")
    print(f"Edges to connect components = k - 1 = {edges_to_connect}")
    print(f"Leaves after connecting = l - 2*(k-1) = {leaves_after_connecting}")
    print(f"Edges for 2-connectivity = ceil({leaves_after_connecting}/2) = {edges_for_2_connectivity}")
    print(f"Total edges = {edges_to_connect} + {edges_for_2_connectivity} = {int(total_edges)}")
    print("\nThe minimal number of new edges to add to G' to make it 2-edge-connected is:")
    print(int(total_edges))

solve()