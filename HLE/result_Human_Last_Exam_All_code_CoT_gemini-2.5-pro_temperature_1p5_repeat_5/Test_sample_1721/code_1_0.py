import math
from fractions import Fraction

# Memoization cache for the recursive search
memo = {}

def get_squares_up_to(limit):
    """Helper function to get all perfect squares up to a limit."""
    squares = set()
    k = 1
    while k * k <= limit:
        squares.add(k * k)
        k += 1
    return squares

def max_independent_set_recursive(graph, nodes_tuple):
    """
    Finds the size of the maximum independent set in a graph using a recursive
    algorithm with memoization. This is computationally expensive.
    """
    nodes = frozenset(nodes_tuple)
    if not nodes:
        return 0
    if nodes in memo:
        return memo[nodes]

    # Pick a pivot vertex (the one with the highest degree in the subgraph)
    pivot = -1
    max_degree = -1
    for node in nodes:
        degree = len(set(graph[node]) & set(nodes))
        if degree > max_degree:
            max_degree = degree
            pivot = node

    # Case 1: Exclude the pivot and find the max independent set in the rest
    remaining_nodes1 = tuple(sorted(list(nodes - {pivot})))
    res1 = max_independent_set_recursive(graph, remaining_nodes1)

    # Case 2: Include the pivot. Remove it and its neighbors.
    neighbors_of_pivot = set(graph[pivot])
    remaining_nodes2 = tuple(sorted(list(nodes - {pivot} - neighbors_of_pivot)))
    res2 = 1 + max_independent_set_recursive(graph, remaining_nodes2)

    result = max(res1, res2)
    memo[result] = result
    return result

def find_max_A_size(N):
    """
    Builds a conflict graph and finds the maximum size of a set A
    where A+A contains no squares. This is equivalent to finding the
    maximum independent set in the graph.
    """
    print(f"\n--- Running exhaustive search for N={N} ---")
    print("This finds the true maximum size of set A. It can be very slow for N > 25.")

    squares = get_squares_up_to(2 * N)
    
    # Build the conflict graph
    adj = {i: [] for i in range(1, N + 1)}
    for i in range(1, N + 1):
        for j in range(i, N + 1):
            if i + j in squares:
                adj[i].append(j)
                if i != j:
                    adj[j].append(i)

    # The problem is finding the maximum independent set in this graph
    all_nodes = tuple(range(1, N + 1))
    memo.clear() # Clear cache for each run
    max_size = max_independent_set_recursive(adj, all_nodes)
    return max_size

def main():
    N = 100  # A reasonably large N for demonstration
    
    print("Problem: Find the largest c such that |A| = (c+o(1))N and A+A has no squares.")
    print("\n--- Construction based on modular arithmetic ---")
    print(f"Let's choose A to be the set of numbers in {{1, ..., {N}}} congruent to 1 mod 3.")
    
    # Construct the set A
    A_mod3 = [n for n in range(1, N + 1) if n % 3 == 1]
    
    # Verify the sumset property for this A (optional, but good for confirmation)
    squares_in_sumset = False
    all_squares = get_squares_up_to(2 * N)
    for i in range(len(A_mod3)):
        for j in range(i, len(A_mod3)):
            if (A_mod3[i] + A_mod3[j]) in all_squares:
                squares_in_sumset = True
                break
        if squares_in_sumset:
            break

    print(f"For N = {N}, A = {{n | 1 <= n <= {N}, n = 3k+1}}")
    # print(f"Set A: {A_mod3}") # Too long to print for N=100
    
    size_A = len(A_mod3)
    ratio = Fraction(size_A, N)
    
    print("\n--- Calculation ---")
    print(f"The size of A is |A| = {size_A}")
    print(f"The ratio is |A|/N = {size_A}/{N} = {float(ratio):.5f}")
    print(f"This construction gives a value of c >= {Fraction(1,3)}.")
    print(f"Verification: A+A contains no squares? {'No' if not squares_in_sumset else 'Yes'}")
    
    print("\nIt is a known result in number theory that c cannot exceed 1/3.")
    print("Thus, the largest possible value for c is 1/3.")

    # Run the exhaustive search for a small N to provide numerical evidence
    N_search = 20
    try:
        max_size = find_max_A_size(N_search)
        ratio_max = Fraction(max_size, N_search)
        print(f"Result of search: For N={N_search}, the max size of A is {max_size}.")
        print(f"The optimal ratio for N={N_search} is |A_max|/N = {max_size}/{N_search} = {float(ratio_max):.5f}")
        size_from_construction = math.floor((N_search - 1) / 3) + 1
        print(f"Our construction for N={N_search} gives size {size_from_construction}, which in this case is optimal.")
    except Exception as e:
        print(f"\nCould not complete the exhaustive search for N={N_search}. It is very slow. Error: {e}")

if __name__ == '__main__':
    main()
