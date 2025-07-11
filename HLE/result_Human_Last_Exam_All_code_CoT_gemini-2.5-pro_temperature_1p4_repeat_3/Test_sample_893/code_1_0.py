def solve():
    """
    This function determines for six classes of partially ordered sets whether they
    always have a maximal element (Y), never have one (N), or if it depends (D).
    """

    # Analysis of each case:
    # A) N: For any H-free graph G, G plus an isolated vertex is a larger H-free graph. No maximal element exists.
    # B) Y: Any non-empty finite subset of real numbers has a maximum, which is a maximal element.
    # C) D: The set of negative integers has a maximal element (-1), but the set of positive integers does not.
    # D) Y: The class of uncountable discrete subsets of R is empty, so the condition is vacuously true.
    # E) Y: The set of sequences with the "subsequence of" relation has maximal elements, e.g., any constant sequence.
    # F) Y: The set of sequences with the "is a subsequence of" relation has a maximum element (a universal sequence), which is maximal.
    
    answer = "NYDYYY"
    print(answer)

solve()