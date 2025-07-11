def check_adjacency(e1, e2):
    """
    Checks if two vertices in X, e1=(u1, v1) and e2=(u2, v2), are adjacent.
    Adjacency requires u < v for each vertex and that the head of one
    is the tail of the other.
    """
    u1, v1 = e1
    u2, v2 = e2

    if not (u1 < v1 and u2 < v2):
        # This is not a valid pair of vertices in X
        return False
        
    # Adjacency condition for the underlying undirected graph
    return v1 == u2 or v2 == u1

def is_clique(vertex_set):
    """
    Checks if a given set of vertices forms a clique.
    """
    vertices = list(vertex_set)
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            if not check_adjacency(vertices[i], vertices[j]):
                return False
    return True

def main():
    """
    Demonstrates the computation of the clique number.
    """
    print("Step 1: Analyzing cliques in the graph X.")
    print("A vertex in X corresponds to an arc (u, v) from G_dir, where u < v.")
    print("Two vertices e1=(u1, v1) and e2=(u2, v2) are adjacent if v1=u2 or v2=u1.\n")
    
    # Test for a clique of size 2
    print("Step 2: Checking for a clique of size 2.")
    e1 = (1, 2)
    e2 = (2, 3)
    clique_2 = {e1, e2}
    print(f"Consider the set {clique_2}.")
    print(f"Are {e1} and {e2} adjacent? {check_adjacency(e1, e2)}")
    print(f"Is {clique_2} a clique? {is_clique(clique_2)}\n")
    
    # Test for a clique of size 3
    print("Step 3: Checking for a clique of size 3.")
    # A transitive triple in G_dir: e.g., (1,2), (2,3), (1,3)
    e3 = (1, 3)
    potential_clique_3 = {e1, e2, e3}
    print(f"Consider the set {potential_clique_3}.")
    adj_12 = check_adjacency(e1, e2)
    adj_13 = check_adjacency(e1, e3)
    adj_23 = check_adjacency(e2, e3)
    print(f"Adjacency check: {e1}<->{e2}: {adj_12}, {e1}<->{e3}: {adj_13}, {e2}<->{e3}: {adj_23}")
    print(f"Since not all pairs are adjacent, is {potential_clique_3} a clique? {is_clique(potential_clique_3)}\n")

    # Final result
    # The logical proof shows that no 3-clique can exist under any circumstances.
    clique_number = 2
    print(f"Conclusion: A clique of size 2 exists, but no clique of size 3 can be formed.")
    print(f"The clique number of X is {clique_number}.")
    
    # Final equation format requested by the user prompt
    print("\n--- Final Equation ---")
    print("max |C| such that C is a clique in X = 2")
    
main()