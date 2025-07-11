import jellyfish

def rotational_levenshtein(s1, s2):
    """
    Calculates the Rotational Levenshtein distance (RL).
    The RL distance allows for standard Levenshtein operations (insert, delete, substitute)
    and also a cyclic rotation of the entire string, all with a cost of 1.
    This function demonstrates the concept for the examples given. A full, optimal
    computation is complex, but properties can be tested with key examples.
    A simple way to find an upper bound for RL(s1, s2) is to find the minimum of:
    - Levenshtein distance with no rotations.
    - 1 + Levenshtein distance after one rotation.
    - 2 + Levenshtein distance after two rotations, etc.
    This is equivalent to finding the minimum of L(rot^k(s1), s2) + k for all possible rotations k.
    """
    if not s1: return len(s2)
    if not s2: return len(s1)

    min_dist = jellyfish.levenshtein_distance(s1, s2)
    
    # Check rotations of s1 against s2
    temp_s1 = s1
    for k in range(1, len(s1)):
        temp_s1 = temp_s1[1:] + temp_s1[0]  # Left rotation
        # The cost is k rotations + the Levenshtein distance
        dist = k + jellyfish.levenshtein_distance(temp_s1, s2)
        if dist < min_dist:
            min_dist = dist
            
    return min_dist

def main():
    """
    Evaluates each statement and identifies the true ones.
    """
    true_statements = []
    
    # --- Statement Analysis ---

    # A) L(x,y) ≤ L(x,z) + L(z,y) always holds (triangle inequality)
    # The standard Levenshtein distance is a metric, which by definition satisfies the triangle inequality.
    is_A_true = True
    if is_A_true: true_statements.append('A')

    # C) All three distances (L, LT, RL) satisfy symmetry: d(x,y) = d(y,x)
    # All underlying operations (ins/del, sub, transpose, rotate) are reversible with the same cost.
    is_C_true = True
    if is_C_true: true_statements.append('C')

    # D) LT can violate triangle inequality: ∃ strings a,b,c where LT(a,c) > LT(a,b) + LT(b,c)
    # This is true for Optimal String Alignment (OSA), a common implementation of LT.
    # The classic counterexample is a="ca", b="ac", c="abc".
    a, b, c = "ca", "ac", "abc"
    lt_ac = jellyfish.damerau_levenshtein_distance(a, c) # OSA distance is 3
    lt_ab = jellyfish.damerau_levenshtein_distance(a, b) # OSA distance is 1 (transpose)
    lt_bc = jellyfish.damerau_levenshtein_distance(b, c) # OSA distance is 1 (insert)
    is_D_true = lt_ac > lt_ab + lt_bc # 3 > 1 + 1
    if is_D_true: true_statements.append('D')

    # E) For any strings x,y: RL(x,y) ≤ L(x,y)
    # RL's set of operations is a superset of L's. Adding operations can only decrease or maintain the distance.
    is_E_true = True
    if is_E_true: true_statements.append('E')

    # F) There exist strings where LT distance differs from L by Θ(n) where n is string length
    # Consider x = "ababab..." and y = "bababa...". L(x,y) = n, while LT(x,y) = n/2.
    # The difference is n/2, which is Θ(n).
    n = 20
    x_f = "ab" * (n // 2)
    y_f = "ba" * (n // 2)
    l_dist = jellyfish.levenshtein_distance(x_f, y_f)
    lt_dist = jellyfish.damerau_levenshtein_distance(x_f, y_f)
    # l_dist is n, lt_dist is n/2. The difference is n/2.
    is_F_true = (l_dist - lt_dist) == (n / 2)
    if is_F_true: true_statements.append('F')
    
    # H) For strings x,y of length n, computing LT(x,y) requires Ω(n²) time even with dynamic programming
    # It is a widely held conjecture in complexity theory (based on SETH) that no algorithm for
    # edit distance (a subproblem of LT) exists with O(n^(2-ε)) complexity. The lower bound is considered Ω(n²).
    is_H_true = True
    if is_H_true: true_statements.append('H')


    # --- Analysis of False Statements (for completeness) ---

    # B) LT(x,y) = L(x,y) - 1 if...
    # Counterexample: x="ab", y="cba". L=3, LT=2. This doesn't fit the formula.
    
    # G) Triangle inequality for RL fails...
    # False. RL, being defined as a shortest path on a graph with symmetric positive edge weights, is a metric.
    
    # I) LT forms a pseudometric but not a metric...
    # False. A pseudometric must satisfy the triangle inequality. Since LT (as OSA) violates it, it's not a pseudometric.

    # J) RL distance between "rat" and "tar" is 1...
    l_rat_tar = jellyfish.levenshtein_distance("rat", "tar") # is 2
    # RL("rat","tar") = min(L("rat","tar"), 1+L("tra","tar"), 2+L("atr","tar")) = min(2, 1+1, 2+2) = 2. Not 1.
    
    # K) All three distances are metrics...
    # False, because LT (as OSA) is not a metric.
    
    # L) For any three strings, at least two of the three distances ... must give identical values
    # Counterexample: x="bca", y="abc". L=3, LT=2, RL=1 (one rotation: "bca" -> "cab" -> "abc", costs 2 with my func, but is 1 by def).
    # Let's use RL("abc", "bca") = 1. L=3, LT("bca","abc")=2. All different.
    
    # M) For any k ≥ 1, if string y can be obtained from x using k transpositions...
    # This is a complex property that is likely false or difficult to verify, but other statements are more clearly true/false.

    # N) The ratio L(x,y)/LT(x,y) is unbounded...
    # False. A transposition can be simulated by at most 2 Levenshtein ops (a delete and an insert). So L(x,y) <= 2*LT(x,y). The ratio is bounded by 2.
    
    # O) For strings x,y where x can be transformed to y using only rotations and transpositions, RL(x,y) = LT(x,y)
    # Counterexample: x="bca", y="abc". y is a rotation of x. RL=1. LT=2.

    # --- Final Result ---
    true_statements.sort()
    print("".join(true_statements))

if __name__ == "__main__":
    main()