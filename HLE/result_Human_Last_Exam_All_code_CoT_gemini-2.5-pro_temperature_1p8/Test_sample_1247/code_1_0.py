import itertools

def count_inversions(p):
    """Counts the number of inversions in a permutation."""
    inv_count = 0
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inv_count += 1
    return inv_count

def has_1324_pattern(p):
    """Checks if a permutation contains the 1324 pattern."""
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    if p[i] < p[k] < p[j] < p[l]:
                        return True
    return False

def calculate_av_n_k(n, k):
    """Calculates the number of permutations of length n with k inversions avoiding 1324."""
    count = 0
    # Generate all permutations of length n
    elements = range(1, n + 1)
    for p in itertools.permutations(elements):
        # Check if the number of inversions is k
        if count_inversions(p) == k:
            # Check if it avoids the 1324 pattern
            if not has_1324_pattern(p):
                count += 1
    return count

# For n=333, direct computation is not feasible.
# The value is constant for n >= 5. We compute for n=5 as a proxy for n=333.
n = 5
k = 3
result = calculate_av_n_k(n, k)
print(result)
