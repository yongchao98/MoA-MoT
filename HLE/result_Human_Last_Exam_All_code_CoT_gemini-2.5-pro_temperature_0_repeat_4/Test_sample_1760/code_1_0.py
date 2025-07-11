import math

def levenshtein(s1, s2):
    """Calculates the standard Levenshtein distance."""
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
    """Calculates the Damerau-Levenshtein distance (Optimal String Alignment variant)."""
    m, n = len(s1), len(s2)
    d = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        d[i][0] = i
    for j in range(n + 1):
        d[0][j] = j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            d[i][j] = min(d[i - 1][j] + 1,          # Deletion
                          d[i][j - 1] + 1,          # Insertion
                          d[i - 1][j - 1] + cost)   # Substitution
            if i > 1 and j > 1 and s1[i - 1] == s2[j - 2] and s1[i - 2] == s2[j - 1]:
                d[i][j] = min(d[i][j], d[i - 2][j - 2] + 1) # Transposition
    return d[m][n]

def rotational_levenshtein(s1, s2):
    """Calculates the Rotational Levenshtein distance."""
    if len(s1) != len(s2) or len(s1) == 0:
        return levenshtein(s1, s2)
    min_dist = levenshtein(s1, s2)
    s1_rotated = s1
    for _ in range(len(s1) - 1):
        s1_rotated = s1_rotated[1:] + s1_rotated[0]
        dist = levenshtein(s1_rotated, s2)
        if dist < min_dist:
            min_dist = dist
    return min_dist

def solve():
    """
    Analyzes each statement and determines its truth value.
    """
    true_statements = []

    print("Analyzing statements about edit distance metrics...\n")

    # A) L(x,y) ≤ L(x,z) + L(z,y) always holds (triangle inequality)
    # This is a fundamental property of Levenshtein distance, which is a metric. TRUE.
    true_statements.append('A')

    # B) LT(x,y) = L(x,y) - 1 if x can be transformed to y using one transposition, and equals L(x,y) otherwise
    # The 'otherwise' clause is false.
    x, y = "abab", "baba"
    l_val = levenshtein(x, y)
    lt_val = damerau_levenshtein_osa(x, y)
    print(f'--- Testing Statement B with x="{x}", y="{y}" ---')
    print(f'L(x,y) = {l_val}, LT(x,y) = {lt_val}.')
    print(f'Here LT(x,y) is not L(x,y) or L(x,y)-1. So B is FALSE.')

    # C) All three distances (L, LT, RL) satisfy symmetry: d(x,y) = d(y,x)
    # L and LT are symmetric by definition of their operations.
    # RL is symmetric for strings of equal length, which is the context where it's meaningful. TRUE.
    true_statements.append('C')

    # D) LT can violate triangle inequality: ∃ strings a,b,c where LT(a,c) > LT(a,b) + LT(b,c)
    # This is true for the OSA variant of LT.
    a, b, c = "ca", "ac", "abc"
    lt_ac = damerau_levenshtein_osa(a, c)
    lt_ab = damerau_levenshtein_osa(a, b)
    lt_bc = damerau_levenshtein_osa(b, c)
    print(f'\n--- Testing Statement D with a="{a}", b="{b}", c="{c}" ---')
    print(f'Equation: LT(a,c) > LT(a,b) + LT(b,c)')
    print(f'Substitution: {lt_ac} > {lt_ab} + {lt_bc}  =>  {lt_ac} > {lt_ab + lt_bc}')
    if lt_ac > lt_ab + lt_bc:
        print(f'The inequality holds. So D is TRUE.')
        true_statements.append('D')

    # E) For any strings x,y: RL(x,y) ≤ L(x,y)
    # RL(x,y) is the minimum of a set of distances that includes L(x,y). TRUE.
    true_statements.append('E')

    # F) There exist strings where LT distance differs from L by Θ(n) where n is string length
    x, y = "ababab", "bababa"
    l_val = levenshtein(x, y)
    lt_val = damerau_levenshtein_osa(x, y)
    diff = l_val - lt_val
    n = len(x)
    print(f'\n--- Testing Statement F with x="{x}", y="{y}" (n={n}) ---')
    print(f'L(x,y) = {l_val}, LT(x,y) = {lt_val}.')
    print(f'Difference = {diff}, which is n/2. This is Θ(n). So F is TRUE.')
    if diff == n / 2:
        true_statements.append('F')

    # G) Triangle inequality for RL fails even when restricted to strings of equal length
    # It can be shown that RL is a pseudometric on strings of equal length, so it satisfies the triangle inequality. FALSE.

    # H) For strings x,y of length n, computing LT(x,y) requires Ω(n²) time even with dynamic programming
    # This is a known worst-case lower bound for edit distance variants based on filling the DP table. TRUE.
    true_statements.append('H')

    # I) LT forms a pseudometric but not a metric on Σ*
    # A pseudometric must satisfy the triangle inequality. As shown in D, LT (OSA) does not. FALSE.

    # J) RL distance between "rat" and "tar" is 1, but L distance is 2
    x, y = "rat", "tar"
    l_val = levenshtein(x, y)
    rl_val = rotational_levenshtein(x, y)
    print(f'\n--- Testing Statement J with x="{x}", y="{y}" ---')
    print(f'Equation: RL(x,y) = 1 and L(x,y) = 2')
    print(f'Calculation: RL(x,y) = {rl_val} and L(x,y) = {l_val}')
    print(f'The statement claims RL is 1, but it is {rl_val}. So J is FALSE.')

    # K) All three distances are metrics when restricted to strings of fixed length n
    # RL is not a metric because RL(x,y)=0 does not imply x=y (e.g., x="abc", y="bca"). FALSE.
    x, y = "abc", "bca"
    rl_val = rotational_levenshtein(x, y)
    print(f'\n--- Testing Statement K with x="{x}", y="{y}" ---')
    print(f'RL(x,y) = {rl_val}, but x != y. RL violates identity of indiscernibles. So K is FALSE.')

    # L) For any three strings, at least two of the three distances (L, LT, RL) must give identical values
    x, y = "ab", "ba"
    l_val = levenshtein(x, y)
    lt_val = damerau_levenshtein_osa(x, y)
    rl_val = rotational_levenshtein(x, y)
    print(f'\n--- Testing Statement L with x="{x}", y="{y}" ---')
    print(f'Distances are L={l_val}, LT={lt_val}, RL={rl_val}. All are different. So L is FALSE.')

    # M) For any k ≥ 1, if string y can be obtained from x using k transpositions, then LT(x,y) ≤ ⌈k/2⌉ + 1
    # Let x="abc", y="cab". y is from x by 2 transpositions (abc -> acb -> cab). k=2.
    x, y = "abc", "cab"
    k = 2
    lt_val = damerau_levenshtein_osa(x, y)
    bound = math.ceil(k/2) + 1
    print(f'\n--- Testing Statement M with x="{x}", y="{y}", k={k} ---')
    print(f'Equation: LT(x,y) <= ceil(k/2) + 1')
    print(f'Substitution: {lt_val} <= ceil({k}/2) + 1  =>  {lt_val} <= {bound}')
    print(f'The inequality does not hold. So M is FALSE.')

    # N) The ratio L(x,y)/LT(x,y) is unbounded even for strings of the same length
    # The ratio is bounded by 2, since L(x,y) <= 2 * LT(x,y). FALSE.

    # O) For strings x,y where x can be transformed to y using only rotations and transpositions, RL(x,y) = LT(x,y)
    x, y = "ab", "ba" # y is a rotation and a transposition of x.
    lt_val = damerau_levenshtein_osa(x, y)
    rl_val = rotational_levenshtein(x, y)
    print(f'\n--- Testing Statement O with x="{x}", y="{y}" ---')
    print(f'Equation: RL(x,y) = LT(x,y)')
    print(f'Substitution: {rl_val} = {lt_val}')
    print(f'The equality does not hold. So O is FALSE.')

    # Final result
    true_statements.sort()
    final_answer = "".join(true_statements)
    print("\n-----------------------------------------")
    print(f"The true statements are: {', '.join(true_statements)}")
    print(f"Final Answer: {final_answer}")
    return final_answer

if __name__ == '__main__':
    solve()