import collections

# Using a dictionary for memoization to store results for subproblems
memo = {}

def count_antichains(s_tuple):
    """
    Recursively counts the number of antichains in a poset defined by divisibility.
    The poset is given by the set of numbers in s_tuple.
    """
    # Base case: The only antichain in an empty set is the empty set itself.
    if not s_tuple:
        return 1
    
    # Use a frozenset for the memoization key as it's hashable and order-independent.
    s_frozenset = frozenset(s_tuple)
    if s_frozenset in memo:
        return memo[s_frozenset]

    s = set(s_tuple)
    
    # Pick an element 'x' to pivot the recursion. 
    # Choosing a maximal element is a good heuristic to simplify subproblems.
    # A maximal element is one that does not divide any other element in the set.
    # We can efficiently find one by picking the largest number in the set.
    x = max(s)

    # --- Recurrence Relation ---
    # a(P) = a(P \ {x}) + a(P \ N[x])
    # where N[x] is the set of all elements comparable to x (including x).

    # Case 1: Antichains that DO NOT contain x.
    # We simply count antichains in the smaller poset without x.
    s_without_x = s.copy()
    s_without_x.remove(x)
    
    # The result for this branch of the recursion
    res = count_antichains(tuple(sorted(list(s_without_x))))

    # Case 2: Antichains that DO contain x.
    # These are of the form {x} U A', where A' is an antichain in the subposet
    # that has no elements comparable to x.
    # N[x] = {y in S | y divides x or x divides y}
    n_x = set()
    for y in s:
        if x % y == 0 or y % x == 0:
            n_x.add(y)
    
    s_prime = s.difference(n_x)
    
    # Add the result from the second branch of the recursion
    res += count_antichains(tuple(sorted(list(s_prime))))
    
    # Memoize the result before returning
    memo[s_frozenset] = res
    return res

# The set S for the problem
S = tuple(range(1, 151))

# The number of open sets in P^-(D_S, tau) is the number of antichains in (D_S, |) plus one.
num_antichains = count_antichains(S)
num_open_sets = num_antichains + 1

# Final equation format as requested
# Number of Open Sets = (Number of Antichains) + 1
print(f"{num_open_sets} = {num_antichains} + 1")