def solve():
    """
    Determines the complexity class for languages G and G'.
    G = {M | exists w forall s: M accepts sw}
    G' = {(k,M) | exists w forall s: |s|=k => M accepts sw}

    Analysis for G:
    The condition is equivalent to "there exists a reset word w that takes every reachable state to the single accept state".
    If such a word w exists, a short one (polynomial in the number of states) is known to exist.
    The problem can be formulated as: Guess a polynomial-length word w (the certificate) and verify it.
    The verification involves:
    1. Finding all reachable states R (in PTIME).
    2. For each state q in R, checking if w takes q to the accept state (in PTIME).
    This structure (existential guess + polynomial-time verifier) is the definition of NP.
    The problem is also NP-hard. Thus, G is in NP.

    Analysis for G':
    The condition is equivalent to "there exists a reset word w that takes every state reachable in exactly k steps to the single accept state".
    Let this set of states be S_k. The problem is resetting the subset S_k.
    The set S_k can be computed in PTIME (using matrix exponentiation).
    However, for an arbitrary-like subset of states, the shortest reset word is not guaranteed to be of polynomial length.
    The problem has an "exists-forall" structure: exists w, forall s of length k.
    This suggests a Sigma_2 complexity.
    The problem can be shown to be Sigma_2-hard by a reduction from QSAT_2, a known Sigma_2-complete problem.
    The membership in Sigma_2 is more complex but holds. It relies on showing that the existence of a (potentially long) reset word can be verified by a non-deterministic machine with an NP oracle.
    Therefore, G' is Sigma_2-complete.
    """
    g_complexity = "NP"
    g_prime_complexity = "Sigma_2"
    
    print(f"{g_complexity}, {g_prime_complexity}")

solve()