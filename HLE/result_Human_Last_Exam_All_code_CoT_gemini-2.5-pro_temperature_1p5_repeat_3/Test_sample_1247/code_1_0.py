import itertools

def count_inversions(perm):
    """Counts the number of inversions in a permutation."""
    inversions = 0
    n = len(perm)
    for i in range(n):
        for j in range(i + 1, n):
            if perm[i] > perm[j]:
                inversions += 1
    return inversions

def has_pattern_1324(perm):
    """Checks if a permutation contains the 1324 pattern."""
    n = len(perm)
    # The pattern 1324 is defined by four indices i < j < k < l
    # such that perm[i] < perm[k] < perm[j] < perm[l].
    # This check is computationally intensive, but feasible for small n.
    for i, j, k, l in itertools.combinations(range(n), 4):
        if perm[i] < perm[k] < perm[j] < perm[l]:
            return True
    return False

def calculate_av_nk_1324(n, k):
    """Calculates av_n^k(1324)."""
    avoiding_perms_count = 0
    # Generate all permutations of {1, 2, ..., n}
    for p in itertools.permutations(range(1, n + 1)):
        # Check if the permutation has k inversions
        if count_inversions(p) == k:
            # Check if the permutation avoids the 1324 pattern
            if not has_pattern_1324(p):
                avoiding_perms_count += 1
    return avoiding_perms_count

# The problem asks for av_333^3(1324).
# As explained, for k=3 and pattern 1324, the value is constant for n >= 5.
# We compute the value for n=5, which will be our answer.
n_compute = 5
k_compute = 3
result = calculate_av_nk_1324(n_compute, k_compute)

# Print the final equation with the requested numbers.
n_final = 333
k_final = 3
pattern_final = "1324"

print(f"av_{n_final}^{k_final}({pattern_final}) = {result}")