import math

def is_affine(s):
    """
    Checks if the relation corresponding to S is affine.
    For a symmetric relation on 4 variables, this corresponds to a known list of 7 sets S.
    """
    affine_sets = [
        set(), {0}, {4}, {0, 4}, {1, 3}, {0, 2, 4}, {0, 1, 2, 3, 4}
    ]
    return s in affine_sets

def is_horn(s):
    """
    Checks if the relation corresponding to S is Horn.
    This translates to a closure property on S:
    For any k, l in S, the integer interval [max(0, k+l-4), min(k,l)] must be a subset of S.
    """
    if not s:
        return True
    for k in s:
        for l in s:
            lower = max(0, k + l - 4)
            upper = min(k, l)
            for m in range(lower, upper + 1):
                if m not in s:
                    return False
    return True

def is_dual_horn(s):
    """
    Checks if the relation corresponding to S is dual-Horn.
    This translates to a closure property on S:
    For any k, l in S, the integer interval [max(k,l), min(4, k+l)] must be a subset of S.
    """
    if not s:
        return True
    for k in s:
        for l in s:
            lower = max(k, l)
            upper = min(4, k + l)
            for m in range(lower, upper + 1):
                if m not in s:
                    return False
    return True

def main():
    """
    Calculates the number of NP-hard cases by subtracting the number of P-time cases from the total.
    """
    p_time_count = 0
    total_count = 2**5
    
    p_time_sets = set()

    for i in range(total_count):
        s = set()
        for j in range(5):
            if (i >> j) & 1:
                s.add(j)

        # The problem is in P if the corresponding relation is affine, Horn, or dual-Horn.
        if is_affine(s) or is_horn(s) or is_dual_horn(s):
            p_time_sets.add(frozenset(s))
    
    p_time_count = len(p_time_sets)
    np_hard_count = total_count - p_time_count

    print(f"Total number of possible sets S: {total_count}")
    print(f"Number of sets S for which the problem is in P: {p_time_count}")
    print(f"Number of sets S for which the problem is NP-hard: {total_count} - {p_time_count} = {np_hard_count}")

if __name__ == '__main__':
    main()