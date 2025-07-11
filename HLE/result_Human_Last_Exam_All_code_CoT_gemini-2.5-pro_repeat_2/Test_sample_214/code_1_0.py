import sys

# It's recommended to increase recursion limit for deep recursive calls.
sys.setrecursionlimit(2000)

# Memoization table to store results for subsets that have already been computed.
memo = {}

def count_antichains(s_tuple):
    """
    Recursively counts the number of antichains in a set of numbers
    under the divisibility relation. An antichain is a subset where no
    element divides another.
    The function uses memoization to store and retrieve results for
    previously seen subsets, avoiding recomputation.
    """
    # Base case: The only antichain in an empty set is the empty set itself.
    if not s_tuple:
        return 1

    # Use a canonical representation (sorted tuple) for the memoization key.
    s_tuple = tuple(sorted(s_tuple))
    if s_tuple in memo:
        return memo[s_tuple]

    # Recursive step: Pick an element 'x' (the largest one) and partition
    # the antichains into two groups: those containing 'x' and those not.
    x = s_tuple[-1]
    rest = s_tuple[:-1]

    # 1. Count antichains that DO NOT contain 'x'.
    # These are simply the antichains of the remaining set.
    count_without_x = count_antichains(rest)

    # 2. Count antichains that DO contain 'x'.
    # These are of the form {x} U C', where C' is an antichain of elements
    # from 'rest' that are not comparable to 'x'. Since all elements in
    # 'rest' are smaller than 'x', we only need to filter out the divisors of 'x'.
    non_divisors = tuple(y for y in rest if x % y != 0)
    count_with_x = count_antichains(non_divisors)

    # The total number of antichains is the sum of the two counts.
    total = count_without_x + count_with_x
    memo[s_tuple] = total
    return total

# The set S' = {2, 3, ..., 150}
S_prime = tuple(range(2, 151))

# Calculate N', the number of antichains in the divisibility poset on S'.
N_prime = count_antichains(S_prime)

# The total number of open sets is N' + 2.
total_open_sets = N_prime + 2

# Print the final result as an equation.
print(f"The number of antichains in the divisibility poset on {{2, ..., 150}} is {N_prime}.")
print(f"The total number of open sets in P^-(D_S, tau) is {N_prime} + 2 = {total_open_sets}")
print(f"<<<{total_open_sets}>>>")