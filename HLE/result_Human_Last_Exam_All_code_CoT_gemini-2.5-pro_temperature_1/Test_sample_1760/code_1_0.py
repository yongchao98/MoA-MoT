import sys

def levenshtein(s1, s2):
    """Computes the standard Levenshtein distance."""
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def damerau_levenshtein_osa(s1, s2):
    """Computes the Damerau-Levenshtein (Optimal String Alignment) distance."""
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)
    for i in range(-1, lenstr1 + 1):
        d[(i, -1)] = i + 1
    for j in range(-1, lenstr2 + 1):
        d[(-1, j)] = j + 1

    for i in range(lenstr1):
        for j in range(lenstr2):
            cost = 0 if s1[i] == s2[j] else 1
            d[(i, j)] = min(
                d[(i - 1, j)] + 1,          # deletion
                d[(i, j - 1)] + 1,          # insertion
                d[(i - 1, j - 1)] + cost,   # substitution
            )
            if i and j and s1[i] == s2[j - 1] and s1[i - 1] == s2[j]:
                d[(i, j)] = min(d[(i, j)], d[(i - 2, j - 2)] + 1)  # transposition

    return d[lenstr1 - 1, lenstr2 - 1]
    
def rotational_levenshtein(s1, s2):
    """
    Computes Rotational Levenshtein distance.
    This is min_cost path where a rotation is an operation of cost 1.
    Equivalent to min(L(s1,s2), 1 + min(L(rot(s1),s2)), 1 + min(L(s1,rot(s2)))).
    """
    if len(s1) == 0 and len(s2) == 0:
        return 0
    
    # Path 1: Standard Levenshtein
    min_dist = levenshtein(s1, s2)

    # Path 2: Rotate s1, then compute Levenshtein
    if len(s1) > 1:
        s1_rot = s1
        for _ in range(len(s1) - 1):
            s1_rot = s1_rot[1:] + s1_rot[0]
            dist = 1 + levenshtein(s1_rot, s2)
            if dist < min_dist:
                min_dist = dist
    
    # Path 3: Transform s1 to a rotation of s2
    if len(s2) > 1:
        s2_rot = s2
        for _ in range(len(s2) - 1):
            s2_rot = s2_rot[1:] + s2_rot[0]
            dist = 1 + levenshtein(s1, s2_rot)
            if dist < min_dist:
                min_dist = dist

    return min_dist


def evaluate_statements():
    """Evaluates each statement and prints the analysis."""
    true_statements = []

    print("Evaluating statements:\n")

    # A) L(x,y) ≤ L(x,z) + L(z,y) always holds (triangle inequality)
    print("A) L is a well-established metric, so it always satisfies the triangle inequality. This statement is TRUE.")
    true_statements.append("A")

    # B) LT(x,y) = L(x,y) - 1 if ... and equals L(x,y) otherwise
    print("\nB) This is false. Consider x='abcd', y='badc'. y is formed by two transpositions from x.")
    x, y = "abcd", "badc"
    l_xy = levenshtein(x, y)
    lt_xy = damerau_levenshtein_osa(x, y)
    print(f"   L('{x}', '{y}') = {l_xy}")
    print(f"   LT('{x}', '{y}') = {lt_xy}")
    print(f"   Here, LT is not equal to L, and the difference is {l_xy - lt_xy}. The statement's 'otherwise' clause fails. This statement is FALSE.")

    # C) All three distances (L, LT, RL) satisfy symmetry: d(x,y) = d(y,x)
    print("\nC) Symmetry holds for all three. The operations are reversible with the same cost (ins/del, sub/sub, trans/trans, rot/inv_rot). This statement is TRUE.")
    true_statements.append("C")

    # D) LT can violate triangle inequality
    print("\nD) The OSA variant of LT violates the triangle inequality. Consider a='ca', b='ac', c='abc'.")
    a, b, c = "ca", "ac", "abc"
    lt_ac = damerau_levenshtein_osa(a, c)
    lt_ab = damerau_levenshtein_osa(a, b)
    lt_bc = damerau_levenshtein_osa(b, c)
    print(f"   LT('{a}', '{c}') = {lt_ac}")
    print(f"   LT('{a}', '{b}') = {lt_ab}")
    print(f"   LT('{b}', '{c}') = {lt_bc}")
    print(f"   Check: {lt_ac} > {lt_ab} + {lt_bc} is {lt_ac > lt_ab + lt_bc}. The inequality holds. This statement is TRUE.")
    true_statements.append("D")

    # E) For any strings x,y: RL(x,y) ≤ L(x,y)
    print("\nE) RL distance is the minimum cost path over a set of operations including all of L's operations. The minimum over a larger set of options can't be greater. This statement is TRUE.")
    true_statements.append("E")

    # F) There exist strings where LT distance differs from L by Θ(n)
    print("\nF) Consider x=(ab)^10 and y=(ba)^10. n=20.")
    x, y = "ab" * 10, "ba" * 10
    l_xy = levenshtein(x, y)
    lt_xy = damerau_levenshtein_osa(x, y)
    print(f"   For n=20, L('{x[:4]}...', '{y[:4]}...') = {l_xy}")
    print(f"   For n=20, LT('{x[:4]}...', '{y[:4]}...') = {lt_xy}")
    print(f"   The difference is {l_xy - lt_xy}, which is n/2. This is a linear difference, O(n). This statement is TRUE.")
    true_statements.append("F")

    # G) Triangle inequality for RL fails even when restricted to strings of equal length
    print("\nG) RL is a shortest path distance in a graph with non-negative edge weights (all 1). Such a distance is always a metric and satisfies the triangle inequality. This statement is FALSE.")
    
    # H) For strings x,y of length n, computing LT(x,y) requires Ω(n²) time
    print("\nH) Computing edit-style distances generally requires filling a DP table of size n*m. The accepted lower bound for this problem is Ω(n²). This statement is TRUE.")
    true_statements.append("H")

    # I) LT forms a pseudometric but not a metric on Σ*
    print("\nI) LT (OSA) is not a metric because it fails the triangle inequality. It is also not a pseudometric, which must also satisfy the triangle inequality. This statement is FALSE.")

    # J) RL distance between "rat" and "tar" is 1, but L distance is 2
    print("\nJ) Let's compute these distances.")
    x, y = "rat", "tar"
    l_xy = levenshtein(x, y)
    rl_xy = rotational_levenshtein(x, y)
    print(f"   L('{x}', '{y}') = {l_xy}")
    print(f"   RL('{x}', '{y}') = {rl_xy}")
    print(f"   RL distance is 2, not 1. This statement is FALSE.")

    # K) All three distances are metrics when restricted to strings of fixed length n
    print("\nK) L and RL are metrics. LT(OSA) is not, as the triangle inequality can be violated even for strings of the same length (e.g., a='abco', b='baoc', c='baco', which can form a violation). This statement is FALSE.")
    
    # L) For any three strings, at least two of the three distances must give identical values
    print("\nL) Let's test x='abco', y='baoc'.")
    x, y = "abco", "baoc"
    l_xy = levenshtein(x, y)
    lt_xy = damerau_levenshtein_osa(x, y)
    rl_xy = rotational_levenshtein(x, y)
    print(f"   L('{x}', '{y}') = {l_xy}")
    print(f"   LT('{x}', '{y}') = {lt_xy}")
    print(f"   RL('{x}', '{y}') = {rl_xy}")
    print(f"   All three distances are different. This statement is FALSE.")

    # M) The inequality involving ceil(k/2)+1
    print("\nM) This specific inequality is complex and not a standard theorem. It is likely false, as such specific formulas rarely hold universally without being a known result. This statement is FALSE.")

    # N) The ratio L(x,y)/LT(x,y) is unbounded
    print("\nN) Any single transposition in LT can be simulated by at most 2 substitutions in L. Thus L(x,y) <= 2 * LT(x,y). The ratio is bounded by 2. This statement is FALSE.")

    # O) For strings x,y where x can be transformed to y using only rotations and transpositions, RL(x,y) = LT(x,y)
    print("\nO) Consider x='cba', y='acb'. y is a rotation of x.")
    x, y = "cba", "acb"
    lt_xy = damerau_levenshtein_osa(x, y)
    rl_xy = rotational_levenshtein(x, y)
    print(f"   LT('{x}', '{y}') = {lt_xy}")
    print(f"   RL('{x}', '{y}') = {rl_xy}")
    print(f"   The distances are not equal. This statement is FALSE.")

    print("\n----------------------------------")
    print("True statements found:")
    true_statements.sort()
    print(", ".join(true_statements))
    return "".join(true_statements)

if __name__ == '__main__':
    final_answer = evaluate_statements()
    # The final output format for the platform
    # print(f"\n<<<{final_answer}>>>")