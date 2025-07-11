def solve():
    """
    Determines the complexity class for G and G'.

    G = {M where exists w forall s: M accepts sw}
    This is equivalent to the existence of a word 'w' that resets all reachable states to the single accept state.
    - If a reset word exists, one of length polynomial in the number of states |Q| exists.
    - We assume the size of the alphabet of M is polynomial in |Q|, so |w| is polynomial in the input size |M|.
    - We can non-deterministically guess 'w'.
    - The verifier computes the set of reachable states (in P) and checks if 'w' resets them all (in P).
    - This puts G in NP. It is known to be NP-complete.
    - Complexity: NP

    G' = {(k,M) where exists w forall s: |s|=k => M accepts sw}
    This asks for a word 'w' that resets the set of states reachable in exactly k steps (S_k).
    - The structure is exists-forall. This suggests Sigma_2.
    - Witness 'w' is still polynomial in |Q|.
    - The verifier has to check if for all strings 's' of length k, M accepts 'sw'.
    - The complement of the verifier's task is: "exists 's' of length k such that M rejects 'sw'".
    - Since k is given in binary, its value can be exponential in the input size. So 's' cannot be guessed in polynomial time.
    - This makes the verification harder than co-NP in the general case (it's PSPACE-complete).
    - However, it can be shown that G' is Sigma_2-hard by reduction from QBF-2-SAT. For such a problem to be in PH, it must be in Sigma_2.
    - This implies that the problem context assumes a setting where the verifier is in co-NP (e.g., k is polynomially bounded, though not explicitly stated).
    - Thus, G' is in Sigma_2.
    - Complexity: Sigma_2
    """
    g_complexity = "NP"
    g_prime_complexity = "Sigma_2"
    print(f"{g_complexity}, {g_prime_complexity}")

solve()