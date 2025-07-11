import math
from itertools import combinations

def solve_count_ans(adj, k):
    """
    Calculates the number of answers for the formula phi_k.

    Args:
        adj (dict): The graph represented as an adjacency list.
        k (int): The integer parameter k for the formula phi_k.
    """
    
    V = list(adj.keys())
    sur_cache = {}

    def combinations_local(n, k_comb):
        if k_comb < 0 or k_comb > n:
            return 0
        if k_comb == 0 or k_comb == n:
            return 1
        if k_comb > n // 2:
            k_comb = n - k_comb
        
        res = 1
        for i in range(k_comb):
            res = res * (n - i) // (i + 1)
        return res

    def sur(k_sur, n_sur):
        """
        Computes the number of surjections from a set of size k to a set of size n.
        """
        if (k_sur, n_sur) in sur_cache:
            return sur_cache[(k_sur, n_sur)]
        
        if k_sur < n_sur:
            return 0
        if n_sur == 0:
            return 1 if k_sur == 0 else 0

        total = 0
        for j in range(n_sur + 1):
            term = combinations_local(n_sur, j) * (j ** k_sur)
            if (n_sur - j) % 2 == 1:
                total -= term
            else:
                total += term
        sur_cache[(k_sur, n_sur)] = total
        return total

    def get_common_neighbors(nodes):
        """
        Finds the set of common neighbors for a given set of nodes.
        """
        if not nodes:
            return set(V)
        
        # Note: using tuple(nodes) for hashability if nodes is a list
        common_set = set(adj.get(nodes[0], []))
        
        for i in range(1, len(nodes)):
            common_set.intersection_update(adj.get(nodes[i], []))
        return common_set

    print(f"Graph vertices: V = {V}")
    print(f"Parameter: k = {k}")
    print("Formula: Sum over subsets B of V (size 1 to k) with common neighbors:")
    print("Total = Sum_{B} sur(k, |B|)")
    print("-" * 30)

    # Reconstructing the sum expression textually for the final output
    sum_parts = []
    final_total = 0
    for b in range(1, k + 1):
        num_sets_with_cn = 0
        for B_tuple in combinations(V, b):
            if get_common_neighbors(list(B_tuple)):
                num_sets_with_cn += 1
        
        if num_sets_with_cn > 0:
            sur_k_b = sur(k, b)
            contribution = num_sets_with_cn * sur_k_b
            sum_parts.append(f"{num_sets_with_cn}*{sur_k_b}")
            final_total += contribution
    
    if not sum_parts:
        print("No answer tuples found.")
        print("Final Answer = 0")
    else:
        equation_str = " + ".join(sum_parts)
        print(f"Total Answers = {equation_str} = {final_total}")
    print("-" * 30)

# Example usage:
# As no graph is given, we use an example graph G.
# G is a path on 4 vertices: 1-2-3-4
example_adj = {
    1: [2],
    2: [1, 3],
    3: [2, 4],
    4: [3]
}
example_k = 3

solve_count_ans(example_adj, example_k)