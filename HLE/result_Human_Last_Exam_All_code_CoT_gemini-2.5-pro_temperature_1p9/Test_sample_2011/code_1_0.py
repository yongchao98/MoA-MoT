import sys

def solve_clique_sizes():
    """
    Determines the maximum possible number of different clique sizes
    that can simultaneously appear as induced subgraphs of a single graph
    on n = 128 vertices.
    """
    
    # The number of vertices in the graph.
    n = 128

    print("Step 1: Understanding the problem.")
    print(f"We have a graph G with n = {n} vertices.")
    print("We want to find the maximum number of different values of 'k' such that an induced K_k (a clique of size k) exists in G.")
    print("-" * 30)

    print("Step 2: Analyzing the set of clique sizes.")
    print("Let ω(G) be the 'clique number' of G, which is the size of the largest clique in G.")
    print("If a graph G has a clique of size ω(G), let's call the set of its vertices 'C'.")
    print("Any subset of C of size 'k' (where 1 <= k <= ω(G)) will also form a clique.")
    print("The subgraph induced by such a k-subset is a K_k.")
    print("This means that if ω(G) is the size of the largest clique, the graph G must contain induced cliques of all sizes: 1, 2, ..., ω(G).")
    print("So, the number of different clique sizes in any graph G is equal to its clique number, ω(G).")
    print("-" * 30)

    print("Step 3: Reframing the problem.")
    print("To maximize the number of different clique sizes, we need to find a graph G that maximizes the clique number ω(G).")
    print("-" * 30)

    print("Step 4: Solving the simplified problem.")
    print(f"For a graph with n = {n} vertices, the size of any clique cannot exceed {n}.")
    print(f"The maximum possible clique number is {n}.")
    print(f"This maximum is achieved by the complete graph, K_{n}, where every vertex is connected to every other vertex.")
    print("-" * 30)
    
    print("Step 5: Final Conclusion and Calculation.")
    # The final equation is simply that the maximum number of sizes equals the number of vertices.
    max_number_of_sizes = n
    
    print(f"The graph K_{n} on n = {n} vertices has a clique number of {n}.")
    print(f"Therefore, it contains induced cliques of sizes 1, 2, ..., {n}.")
    print(f"The total number of different clique sizes is {n}.")
    print(f"\nFinal Equation: max_sizes = n")
    print(f"Final Answer: max_sizes = {max_number_of_sizes}")

if __name__ == "__main__":
    solve_clique_sizes()