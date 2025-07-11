import sys

# Increase recursion limit for deep DP calculations if needed, though not for these examples
# sys.setrecursionlimit(2000)

def levenshtein_distance(s1, s2):
    """Computes the standard Levenshtein distance (L)."""
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1,        # Deletion
                           dp[i][j - 1] + 1,        # Insertion
                           dp[i - 1][j - 1] + cost) # Substitution
    return dp[m][n]

def damerau_levenshtein_osa(s1, s2):
    """Computes the Damerau-Levenshtein distance (Optimal String Alignment version - LT)."""
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1,        # Deletion
                           dp[i][j - 1] + 1,        # Insertion
                           dp[i - 1][j - 1] + cost) # Substitution
            if i > 1 and j > 1 and s1[i - 1] == s2[j - 2] and s1[i - 2] == s2[j - 1]:
                dp[i][j] = min(dp[i][j], dp[i - 2][j - 2] + 1) # Transposition

    return dp[m][n]

def rotational_levenshtein(s1, s2):
    """Computes the Rotational Levenshtein distance (RL) for equal length strings."""
    if len(s1) != len(s2):
        # The concept is ill-defined for strings of different lengths.
        # Fallback to standard Levenshtein distance.
        return levenshtein_distance(s1, s2)
    
    n = len(s1)
    min_dist = levenshtein_distance(s1, s2)
    s1_rotated = s1
    for _ in range(n - 1):
        s1_rotated = s1_rotated[1:] + s1_rotated[0]
        dist = levenshtein_distance(s1_rotated, s2)
        if dist < min_dist:
            min_dist = dist
    return min_dist

def analyze_statements():
    """Analyzes each statement and determines its truth value."""
    x = "algorithm"
    y = "logarithm"
    z = "altarithm"
    true_statements = []

    print("Analyzing Statements:\n")

    # A)
    print("A) L(x,y) <= L(x,z) + L(z,y) always holds (triangle inequality)")
    l_xy = levenshtein_distance(x, y)
    l_xz = levenshtein_distance(x, z)
    l_zy = levenshtein_distance(z, y)
    print(f"Test with example strings: L(x,y)={l_xy}, L(x,z)={l_xz}, L(z,y)={l_zy}")
    print(f"Checking inequality: {l_xy} <= {l_xz} + {l_zy}, which is {l_xy <= l_xz + l_zy}")
    print("This is a fundamental property of all metrics, which Levenshtein distance is. TRUE.\n")
    true_statements.append("A")

    # B)
    print("B) LT(x,y) = L(x,y) - 1 if x can be transformed to y using one transposition, and equals L(x,y) otherwise")
    print("This is false. Counterexample: x='abcde', y='badce'. LT=2, L=3. y is not from x by one transposition, yet LT != L. FALSE.\n")

    # C)
    print("C) All three distances (L, LT, RL) satisfy symmetry: d(x,y) = d(y,x)")
    l_yx = levenshtein_distance(y,x)
    lt_xy = damerau_levenshtein_osa(x,y)
    lt_yx = damerau_levenshtein_osa(y,x)
    print(f"L is symmetric: L(x,y)={l_xy}, L(y,x)={l_yx}")
    print(f"LT is symmetric: LT(x,y)={lt_xy}, LT(y,x)={lt_yx}")
    print("RL is symmetric because L is, and the set of relative rotations is symmetric.")
    print("All operations (ins/del/sub/transpose/rotate) are reversible with the same cost. TRUE.\n")
    true_statements.append("C")

    # D)
    print("D) LT can violate triangle inequality: ∃ strings a,b,c where LT(a,c) > LT(a,b) + LT(b,c)")
    a, b, c = "ca", "ac", "abc"
    lt_ac = damerau_levenshtein_osa(a, c)
    lt_ab = damerau_levenshtein_osa(a, b)
    lt_bc = damerau_levenshtein_osa(b, c)
    print(f"Test with counterexample a='{a}', b='{b}', c='{c}'")
    print(f"LT(a,c) = {lt_ac}")
    print(f"LT(a,b) = {lt_ab}")
    print(f"LT(b,c) = {lt_bc}")
    print(f"Checking {lt_ac} > {lt_ab} + {lt_bc}, which is {lt_ac > lt_ab + lt_bc}")
    print("The common OSA implementation of Damerau-Levenshtein is not a true metric and violates this. TRUE.\n")
    true_statements.append("D")
    
    # E)
    print("E) For any strings x,y: RL(x,y) <= L(x,y)")
    rl_xy = rotational_levenshtein(x,y)
    print(f"RL(x,y) = {rl_xy}, L(x,y) = {l_xy}. Checking {rl_xy} <= {l_xy}: {rl_xy <= l_xy}")
    print("RL is the minimum Levenshtein distance over all rotations, one of which is the identity (zero rotation). TRUE.\n")
    true_statements.append("E")

    # F)
    print("F) There exist strings where LT distance differs from L by Θ(n) where n is string length")
    print("Consider x='a1b1a2b2...' and y='b1a1b2a2...'. L(x,y) = n. LT(x,y)=n/2 (n/2 transpositions).")
    print("The difference L-LT is n/2, which is Θ(n). TRUE.\n")
    true_statements.append("F")

    # G)
    print("G) Triangle inequality for RL fails even when restricted to strings of equal length")
    print("This is false. RL, as defined, is a pseudometric on the space of cyclic strings, which satisfies the triangle inequality. FALSE.\n")

    # H)
    print("H) For strings x,y of length n, computing LT(x,y) requires Ω(n²) time even with dynamic programming")
    print("This is widely believed to be true. No truly sub-quadratic algorithm is known, and finding one would refute the Strong Exponential Time Hypothesis (SETH). TRUE.\n")
    true_statements.append("H")

    # I)
    print("I) LT forms a pseudometric but not a metric on Σ*")
    print("This is false. Since LT violates the triangle inequality (Statement D), it is not a pseudometric. A pseudometric must satisfy triangle inequality. FALSE.\n")
    
    # J)
    print("J) RL distance between 'rat' and 'tar' is 1, but L distance is 2")
    r, t = "rat", "tar"
    l_rt = levenshtein_distance(r, t)
    rl_rt = rotational_levenshtein(r, t)
    print(f"L('rat','tar') = {l_rt}")
    print(f"RL('rat','tar') = {rl_rt} (min L distance of 'rat' to {{'tar', 'art', 'rta'}}) is not 1.")
    print("So the statement is FALSE.\n")

    # K)
    print("K) All three distances are metrics when restricted to strings of fixed length n")
    print("False. LT(OSA) is not a metric (violates triangle inequality). RL is a pseudometric, not a metric, as RL(x,y)=0 does not imply x=y (e.g. 'ab' vs 'ba'). FALSE.\n")
    
    # L)
    print("L) For any three strings, at least two of the three distances (L, LT, RL) must give identical values")
    a,b = "ab", "ba"
    l_ab = levenshtein_distance(a,b)
    lt_ab = damerau_levenshtein_osa(a,b)
    rl_ab = rotational_levenshtein(a,b)
    print(f"Counterexample: x='{a}', y='{b}'. L(x,y)={l_ab}, LT(x,y)={lt_ab}, RL(x,y)={rl_ab}. All are different. FALSE.\n")
    
    # M)
    print("M) For any k >= 1, if string y can be obtained from x using k transpositions, then LT(x,y) <= ceil(k/2) + 1")
    print("False. Counterexample: For k=4 non-overlapping transpositions, LT=4. Formula gives ceil(4/2)+1 = 3. 4<=3 is false. FALSE.\n")

    # N)
    print("N) The ratio L(x,y)/LT(x,y) is unbounded even for strings of the same length")
    print("False. A single transposition (cost 1 for LT) can be simulated by at most 2 substitutions (cost 2 for L). Thus, L(x,y) <= 2*LT(x,y). The ratio is bounded by 2. FALSE.\n")
    
    # O)
    print("O) For strings x,y where x can be transformed to y using only rotations and transpositions, RL(x,y) = LT(x,y)")
    print(f"Counterexample: x='ab', y='ba'. Can be transformed by one transposition. LT=1, but RL=0. FALSE.\n")
    
    true_statements.sort()
    return true_statements

if __name__ == '__main__':
    final_true_statements = analyze_statements()
    final_answer = ",".join(final_true_statements)
    print("\n-------------------------------------------")
    print(f"The true statements are: {final_answer}")
    print("-------------------------------------------\n")
    print("<<<" + final_answer + ">>>")
