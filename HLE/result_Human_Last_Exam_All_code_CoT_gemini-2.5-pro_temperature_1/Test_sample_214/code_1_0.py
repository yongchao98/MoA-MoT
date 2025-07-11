import sys

# It's better to increase recursion limit for deep recursive calls.
sys.setrecursionlimit(2000)

memo = {}

def count_antichains(nodes):
    """
    Recursively counts the number of antichains in a poset defined by divisibility.
    The set of nodes is represented by a sorted tuple to enable memoization.
    """
    # Base case: The empty set has one antichain, the empty antichain.
    if not nodes:
        return 1
    
    nodes_tuple = tuple(nodes)
    # Return stored result if already computed
    if nodes_tuple in memo:
        return memo[nodes_tuple]

    # Pick the largest element 'x' for the recursive step.
    # This is an optimization heuristic.
    x = nodes[-1]
    
    # Group 1: Antichains that do not contain x.
    # These are all antichains of the sub-poset without x.
    nodes_without_x = nodes[:-1]
    count1 = count_antichains(nodes_without_x)

    # Group 2: Antichains that contain x.
    # These are of the form {x} U A, where A is an antichain of elements
    # from the remaining set that are incomparable to x.
    incomparable_nodes = []
    for y in nodes_without_x:
        # Since y < x, y can't be a multiple of x.
        # We only need to check if y divides x.
        if x % y != 0:
            incomparable_nodes.append(y)
    
    count2 = count_antichains(incomparable_nodes)

    # Total number of antichains is the sum of the two groups.
    total = count1 + count2
    memo[nodes_tuple] = total
    
    return total

if __name__ == '__main__':
    # The set S = {1, 2, ..., 150}
    s_nodes = list(range(1, 151))
    
    # The number of open sets is the number of antichains in the divisibility poset on S.
    num_open_sets = count_antichains(s_nodes)
    
    print("The number of open sets in (D_S, tau) is the number of antichains in the divisibility poset on S.")
    print(f"The number of antichains in the divisibility poset on S = {{1, 2, ..., 150}} is: {num_open_sets}")
    print("Based on the derivation, this is also the number of open sets in P^-(D_S, tau).")
    # The final equation is simply Result = num_open_sets
    print(f"Final Answer = {num_open_sets}")