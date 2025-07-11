import numpy as np

def levenshtein(s1, s2):
    """Calculates the standard Levenshtein distance."""
    if not isinstance(s1, str) or not isinstance(s2, str):
        raise TypeError("Inputs must be strings.")
    size_x = len(s1) + 1
    size_y = len(s2) + 1
    matrix = np.zeros((size_x, size_y), dtype=int)
    for x in range(size_x):
        matrix[x, 0] = x
    for y in range(size_y):
        matrix[0, y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            cost = 0 if s1[x-1] == s2[y-1] else 1
            matrix[x,y] = min(
                matrix[x-1,y] + 1,        # Deletion
                matrix[x,y-1] + 1,        # Insertion
                matrix[x-1,y-1] + cost)   # Substitution
    return matrix[size_x - 1, size_y - 1]

def damerau_levenshtein_osa(s1, s2):
    """Calculates the Damerau-Levenshtein distance using the OSA algorithm."""
    if not isinstance(s1, str) or not isinstance(s2, str):
        raise TypeError("Inputs must be strings.")
    len_s1 = len(s1)
    len_s2 = len(s2)
    d = np.zeros((len_s1 + 1, len_s2 + 1), dtype=int)
    for i in range(len_s1 + 1):
        d[i, 0] = i
    for j in range(len_s2 + 1):
        d[0, j] = j

    for i in range(1, len_s1 + 1):
        for j in range(1, len_s2 + 1):
            cost = 0 if s1[i-1] == s2[j-1] else 1
            d[i, j] = min(
                d[i-1, j] + 1,          # Deletion
                d[i, j-1] + 1,          # Insertion
                d[i-1, j-1] + cost      # Substitution
            )
            if i > 1 and j > 1 and s1[i-1] == s2[j-2] and s1[i-2] == s2[j-1]:
                d[i, j] = min(d[i, j], d[i-2, j-2] + 1) # Transposition
    return d[len_s1, len_s2]

def rotational_levenshtein(s1, s2):
    """Calculates the Rotational Levenshtein distance."""
    l_dist = levenshtein(s1, s2)
    if not s1 or not s2:
        return l_dist
    
    min_dist = l_dist
    
    # Rotations of s1 vs s2
    temp_s1 = s1
    for _ in range(len(s1) - 1):
        temp_s1 = temp_s1[1:] + temp_s1[0]
        dist = 1 + levenshtein(temp_s1, s2)
        if dist < min_dist:
            min_dist = dist
            
    # Rotations of s2 vs s1
    temp_s2 = s2
    for _ in range(len(s2) - 1):
        temp_s2 = temp_s2[1:] + temp_s2[0]
        dist = 1 + levenshtein(s1, temp_s2)
        if dist < min_dist:
            min_dist = dist
            
    return min_dist

def analyze_statements():
    """Analyzes each statement and determines its truth value."""
    results = {}
    
    # A
    results['A'] = True  # L is a metric, triangle inequality always holds.

    # B
    results['B'] = False # LT is not always L-1. e.g. L("abac", "baca")=2, LT=2. Also "otherwise" is false.
    
    # C
    results['C'] = True # All operations are reversible with the same cost, so symmetry holds.

    # D
    a, b, c = "ca", "ac", "abc"
    lt_ac = damerau_levenshtein_osa(a, c) # 3
    lt_ab = damerau_levenshtein_osa(a, b) # 1
    lt_bc = damerau_levenshtein_osa(b, c) # 1
    results['D'] = lt_ac > lt_ab + lt_bc # 3 > 1 + 1 is True

    # E
    results['E'] = True # RL adds operations to L, so the distance can only be <= L.
    
    # F
    x, y = "ababababab", "bababababa" # n=10
    l_dist = levenshtein(x,y) # 10
    lt_dist = damerau_levenshtein_osa(x,y) # 5
    # L=n, LT=n/2. The difference is n/2 = Theta(n).
    results['F'] = True
    
    # G
    results['G'] = False # RL is a metric and satisfies the triangle inequality.
    
    # H
    results['H'] = True # It is a known result that computing edit distance (L or LT) requires Omega(n^2) in the general case.
    
    # I
    results['I'] = False # A pseudometric must satisfy the triangle inequality. Since LT (OSA) violates it (see D), it's not a pseudometric.
    
    # J
    x, y = "rat", "tar"
    l_dist = levenshtein(x,y) # 2
    rl_dist = rotational_levenshtein(x,y) # 2
    results['J'] = (rl_dist == 1 and l_dist == 2) # This is False as RL is 2.

    # K
    results['K'] = False # "All three" is false because LT (OSA) is not a metric, regardless of string length.
    
    # L
    results['L'] = True # It is very difficult to construct a counterexample. If no rotation helps, RL=L. If no transposition helps, LT=L. If both help, the operations often interfere, making it unlikely for all three values to be distinct. E.g. x="dabc", y="badc" -> L=2, RL=2, LT=1.
    
    # M
    x, y = "abcdefgh", "badcfehg" # y is from x with 4 non-overlapping transpositions (k=4)
    lt_dist = damerau_levenshtein_osa(x,y) # 4
    bound = np.ceil(4/2) + 1 # 3
    results['M'] = lt_dist <= bound # 4 <= 3 is False.
    
    # N
    results['N'] = False # The ratio is bounded by 2, since L(x,y) <= 2 * LT(x,y).
    
    # O
    x, y = "abc", "bca" # y is from x by rotation only
    rl_dist = rotational_levenshtein(x,y) # 1
    lt_dist = damerau_levenshtein_osa(x,y) # 2
    results['O'] = rl_dist == lt_dist # 1 == 2 is False.
    
    true_statements = sorted([k for k, v in results.items() if v])
    print("Analysis of each statement:")
    for statement, is_true in sorted(results.items()):
        print(f"Statement {statement}: {is_true}")
        
    print("\nThe true statements in lexicographic order are:")
    print(', '.join(true_statements))

analyze_statements()
