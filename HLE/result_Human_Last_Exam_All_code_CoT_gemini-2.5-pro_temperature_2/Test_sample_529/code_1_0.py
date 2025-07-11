import collections
from itertools import combinations
from math import factorial

# Memoization for Stirling numbers of the second kind
_stirling_cache = {}
def stirling_second(n, k):
    """
    Computes the Stirling number of the second kind S2(n, k).
    This is the number of ways to partition a set of n elements into k non-empty subsets.
    """
    if (n, k) in _stirling_cache:
        return _stirling_cache[(n, k)]
    if n == k == 0:
        return 1
    if k == 0 or n == 0 or k > n:
        return 0
    if k == 1 or n == k:
        return 1
    
    result = k * stirling_second(n - 1, k) + stirling_second(n - 1, k - 1)
    _stirling_cache[(n, k)] = result
    return result

def solve_count_ans(graph, k):
    """
    Solves the CountAns problem for a given graph and parameter k.

    Args:
        graph (dict): Adjacency list representation of the graph.
        k (int): The parameter k for the formula phi_k.

    Returns:
        int: The number of answers of phi_k in G.
    """
    nodes = list(graph.keys())
    n = len(nodes)
    total_answers = 0
    
    print("Graph G has nodes:", nodes)
    print(f"Parameter k = {k}\n")
    print("Formula: CountAns(G, k) = sum_{j=1 to k} count_j * j! * S2(k, j)")
    print("-" * 50)
    
    # max_img_size is k because a k-tuple has at most k distinct vertices.
    # A tuple (v1,...,vk) where all v_i are the same has image size 1.
    for j in range(1, k + 1):
        if j > n:
            break
            
        count_j = 0 # Number of subsets of size j with a common neighbor
        
        # Iterate over all subsets of V(G) of size j
        for s_tuple in combinations(nodes, j):
            s_set = set(s_tuple)
            
            # Find the common neighbors of the set S
            if not s_set:
                continue

            # Start with the neighbors of the first vertex in S
            first_vertex = s_tuple[0]
            common_neighbors = set(graph[first_vertex])

            # Intersect with neighbors of other vertices in S
            for i in range(1, len(s_tuple)):
                vertex = s_tuple[i]
                common_neighbors.intersection_update(graph[vertex])
            
            # If there is at least one common neighbor
            if common_neighbors:
                count_j += 1
        
        if count_j > 0:
            fact_j = factorial(j)
            s2_k_j = stirling_second(k, j)
            term = count_j * fact_j * s2_k_j
            
            print(f"For j = {j}:")
            print(f"  - Number of {j}-sets with common neighbors (count_{j}): {count_j}")
            print(f"  - j! = {fact_j}")
            print(f"  - Stirling number S2({k},{j}) = {s2_k_j}")
            print(f"  - Term for j={j}: {count_j} * {fact_j} * {s2_k_j} = {term}")
            print("-" * 50)
            
            total_answers += term

    print(f"\nTotal number of answers = {total_answers}")
    return total_answers

if __name__ == '__main__':
    # Example Usage: A cycle graph C4
    # V = {0, 1, 2, 3}
    # E = {(0,1), (1,2), (2,3), (3,0)}
    # Let's check with k=2
    # The pairs with a common neighbor are:
    # {0,2} have common neighbors {1,3}
    # {1,3} have common neighbors {0,2}
    # So valid tuples are those formed from {0,2} and {1,3}.
    # Tuples from {0,2}: (0,0), (0,2), (2,0), (2,2)
    # Tuples from {1,3}: (1,1), (1,3), (3,1), (3,3)
    # Also tuples from {0}, {1}, {2}, {3} alone work.
    # E.g. (0,0) is an answer because N(0) is not empty. Wait, that's not the definition.
    # (0,0) needs CN({0}) = N(0) to be non-empty. Which it is.
    # Total answers should be sum of terms.
    # j=1: sets {0},{1},{2},{3}. all have non-empty N(v). count_1=4. Term=4*1!*S2(2,1)=4.
    # j=2: sets {0,2}, {1,3}. count_2=2. Term=2*2!*S2(2,2)=4.
    # Total = 4 + 4 = 8.
    
    # A graph can be represented as an adjacency list (dictionary in Python)
    G = {
        0: [1, 3],
        1: [0, 2],
        2: [1, 3],
        3: [0, 2]
    }
    k = 2
    
    solve_count_ans(G, k)
